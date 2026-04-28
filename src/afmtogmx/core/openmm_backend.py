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
            Configuration keys to update. Available keys:

            - ``incl_mol`` (list): molecule names to include (default: all).
            - ``molname_translations`` (dict): maps molecule names from the
              ``.off`` file to PDB-compatible residue names, e.g.
              ``{'H2OQM': 'SOL', 'CYCQM': 'CYC'}``.

        Returns
        -------
        OpenMMBackend
            Returns ``self`` for method chaining.

        Examples
        --------
        >>> off.openmm.set_config(molname_translations={'H2OQM': 'SOL'})
        >>> off.openmm.gen_xml(output='forcefield.xml')
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

        Examples
        --------
        >>> off.openmm.get_config('incl_mol')
        []
        >>> off.openmm.get_config()
        {'incl_mol': [], 'molname_translations': {}}
        """
        if key is None:
            return self.config.copy()
        return self.config.get(key)

    def gen_xml(self, output='forcefield.xml', incl_mol=None, molname_translations=None):
        """Generate an OpenMM XML force field file from parsed ``.off`` data.

        Reads nonbonded, bonded, and charge data from the parent
        :class:`ReadOFF` instance and writes a complete OpenMM-compatible
        ``<ForceField>`` XML file.

        Parameters
        ----------
        output : str, optional
            Path for the output ``.xml`` file. Defaults to
            ``'forcefield.xml'``.
        incl_mol : list of str, optional
            Molecule names to include. If ``None`` or empty, all molecules
            are included. Defaults to ``[]`` (from ``self.config``).
        molname_translations : dict, optional
            Maps molecule names from the ``.off`` file to PDB-compatible
            3-letter residue names, e.g. ``{'H2OQM': 'SOL'}``. Atoms and
            residues in the XML will use the translated names. Defaults to
            ``{}`` (from ``self.config``).

        Returns
        -------
        None

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

        raise NotImplementedError(
            "OpenMMBackend.gen_xml() is not yet implemented. "
            "Use offtoopenmm/off_to_openmm.py in the meantime."
        )
