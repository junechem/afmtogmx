from afmtogmx.core import xml_generation
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
