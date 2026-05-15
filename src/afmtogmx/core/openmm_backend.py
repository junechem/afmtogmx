from afmtogmx.core import xml_generation, pdb_processing
from types import SimpleNamespace


class OpenMMBackend:
    """OpenMM output backend for ReadOFF.

    Accessed via ``off.openmm`` on a :class:`ReadOFF` instance. Owns all
    OpenMM-specific configuration and output generation. Reads shared
    force-field data from the parent via ``self._parent.nonbonded``,
    ``self._parent.bonded``, and ``self._parent.charges``.
    """

    def __init__(self, parent):
        self._parent = parent
        self.config = {
            'incl_mol': [],
            'molname_translations': {},
        }

    def set_config(self, **kwargs):
        """Update the internal configuration with new default values.

        Parameters
        ----------
        **kwargs
            Configuration keys to update.

        Returns
        -------
        OpenMMBackend
            ``self`` for method chaining.

        Notes
        -----
        Available configuration keys and their default values:
        - ``incl_mol`` (list): ``[]``  — molecule names to include (empty = all)
        - ``molname_translations`` (dict): ``{}``  — maps ``.off`` molecule
          names to PDB-compatible residue names, e.g. ``{'H2OQM': 'SOL'}``.

        Examples
        --------
        >>> off = ReadOFF('forcefield.off')
        >>> off.openmm.set_config(molname_translations={'H2OQM': 'SOL'})
        >>> off.openmm.gen_xml()
        """
        self.config.update(kwargs)
        return self

    def get_config(self, key=None):
        """Retrieve configuration value(s).

        Parameters
        ----------
        key : str, optional
            Configuration key to retrieve. If ``None``, returns a copy of
            the full configuration dictionary.

        Returns
        -------
        value or dict
        """
        if key is None:
            return self.config.copy()
        return self.config.get(key)

    def gen_xml(self, output='forcefield.xml', incl_mol=None, molname_translations=None):
        """Generate an OpenMM XML force-field file from parsed ``.off`` data.

        Reads bonded, nonbonded, and charge data from the parent
        :class:`ReadOFF` and writes a complete ``<ForceField>`` XML.
        Atom types are namespaced as ``"<MOLNAME>_<TYPE>"`` so that two
        ``.off`` molecules reusing the same raw type label do not collide.

        Parameters
        ----------
        output : str, optional
            Path for the output ``.xml`` file. Defaults to ``'forcefield.xml'``.
        incl_mol : list of str, optional
            Molecule names to include. Empty/``None`` includes every molecule
            in the parent. Defaults to ``[]`` (from ``self.config``).
        molname_translations : dict, optional
            Maps ``.off`` molecule names to PDB-compatible residue names,
            e.g. ``{'H2OQM': 'SOL'}``. Defaults to ``{}`` (from ``self.config``).

        Returns
        -------
        None

        Notes
        -----
        - Supported nonbonded interaction types: EXP, STR/STRC, SRD, POW, BUC.
          POW is folded into the SRD force (r0=0). BUC is split into EXP+SRD.
        - Unit conversions: kcal/mol → kJ/mol (×4.184), Å → nm (×0.1).
        - Charges live on ``self._parent.charges`` keyed by atom *name*; the
          builder maps name → type via the ATO section before writing the
          ``<NonbondedForce>``.

        Examples
        --------
        >>> off = ReadOFF('forcefield.off')
        >>> off.load_charges_from_file('charges.txt')
        >>> off.openmm.set_config(molname_translations={'H2OQM': 'SOL'})
        >>> off.openmm.gen_xml(output='forcefield.xml')
        """
        # Resolve parameters: explicit value → config → default
        p = SimpleNamespace(**{
            k: v if v is not None else self.config[k]
            for k, v in locals().items() if k in self.config
        })

        bonded    = self._parent.bonded
        nonbonded = self._parent.nonbonded
        charges   = self._parent.charges

        mol_names      = [m for m in bonded if not p.incl_mol or m in p.incl_mol]
        atom_types     = xml_generation.collect_atom_types(bonded, mol_names)
        type_to_charge = xml_generation.build_type_to_charge(bonded, charges, atom_types)

        sections = []
        sections.append(xml_generation.gen_atomtypes(bonded, atom_types))
        sections.append(xml_generation.gen_residues(bonded, mol_names, p.molname_translations))
        sections.append(xml_generation.gen_nonbonded_force(atom_types, type_to_charge))

        for builder in (xml_generation.gen_bond_force,
                        xml_generation.gen_angle_force,
                        xml_generation.gen_dihedral_force):
            section = builder(bonded, mol_names)
            if section:
                sections.append(section)

        exp_entries, str_entries, srd_by_power = xml_generation.collect_nonbonded(
            nonbonded, atom_types,
        )
        if exp_entries:
            sections.append(xml_generation.gen_exp_force(exp_entries, atom_types))
        for power in sorted(srd_by_power.keys()):
            sections.append(xml_generation.gen_srd_force(srd_by_power[power], power, atom_types))
        if str_entries:
            sections.append(xml_generation.gen_str_force(str_entries, atom_types))

        xml_generation.write_xml(output, sections)

    def preprocess_pdb(self, input_pdb, xml_file, output_pdb='output.pdb',
                       maxwarn=0):
        """Preprocess a PDB file for OpenMM compatibility using a force-field XML.

        Renames PDB atoms to match the atom names in the XML's ``<Residues>``
        section and emits CONECT records derived from the XML's bond
        definitions. The XML file is the only source of force-field
        information used; no in-memory ``ReadOFF`` state is consulted, so a
        single XML can be reused to process many PDB files.

        Atom identification uses one of two strategies, chosen automatically:

        - **Topology matching** — used when the input PDB has CONECT records.
          Each residue's bond graph is matched to the XML residue graph by
          element-labelled graph isomorphism, so atoms are identified by
          *connectivity* regardless of their order in the file. Atoms are
          then reordered into the canonical XML order. This is the robust
          path for PDBs whose atom order differs from the ``.off``/XML (e.g.
          PDBs built from a CIF).
        - **Positional matching** — fallback used when the PDB has no CONECT
          records. The Nth atom of each residue receives the Nth XML name;
          the ``maxwarn`` failsafe (below) guards against mis-ordered input.

        Parameters
        ----------
        input_pdb : str
            Path to input PDB file.
        xml_file : str
            Path to an OpenMM-style force-field XML, e.g. one produced by
            :meth:`gen_xml`.
        output_pdb : str, optional
            Path for output PDB file. Defaults to ``'output.pdb'``.
        maxwarn : int, optional
            Number of element-mismatch warnings tolerated on the positional
            (no-CONECT) path before the run aborts. Default ``0`` (strict).
            Ignored when the PDB has CONECT records.

        Returns
        -------
        OpenMMBackend
            ``self`` for method chaining.

        Raises
        ------
        FileNotFoundError
            If ``input_pdb`` or ``xml_file`` does not exist.
        IOError
            If ``output_pdb`` cannot be written.
        SystemExit
            If a residue cannot be matched to the XML by bond topology, or
            the positional element-mismatch count exceeds ``maxwarn``.

        Notes
        -----
        - All existing CONECT records in the input PDB are replaced with
          new ones derived from the XML bond definitions, and atom serial
          numbers are renumbered sequentially.
        - On the topology path, residue atoms are reordered into the XML
          order; on the positional path the original order is kept.
        - Virtual sites (mass=0) are included in the atom renaming.

        Examples
        --------
        >>> off = ReadOFF('butanol_water.off')
        >>> off.change_molecule('H2OQM', 'BLYPSP-4F', ref_mol_name='H2OQM')
        >>> off.openmm.set_config(incl_mol=['H2OQM'], molname_translations={'H2OQM': 'SOL'})
        >>> off.openmm.gen_xml(output='forcefield.xml')
        >>> off.openmm.preprocess_pdb('conf.pdb', 'forcefield.xml', 'conf_processed.pdb')

        Standalone usage (no ReadOFF object needed)::

            from afmtogmx.core import pdb_processing
            topology = pdb_processing.build_residue_topology_from_xml('forcefield.xml')
            pdb_processing.process_pdb('conf.pdb', 'conf_processed.pdb', topology)
        """
        residue_topology = pdb_processing.build_residue_topology_from_xml(xml_file)
        pdb_processing.process_pdb(input_pdb, output_pdb, residue_topology,
                                   maxwarn=maxwarn)
        return self
