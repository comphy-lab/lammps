AMBER to LAMMPS Tutorial
========================

**written by Arun Srikanth Sridhar**

`AMBER2LAMMPS <amber2lammps_>`_ is a Python utility that converts AMBER topology (``.prmtop``),
force-field (``.frcmod``), and coordinate/charge (``.mol2``) files into
LAMMPS data and parameter files. It provides both a command-line interface
and a Python API with built-in validation. This document provides a detailed workflow with examples using the
`AMBER2LAMMPS <amber2lammps_>`_ python utility. The legacy scripts
previously distributed in ``tools/amber2lmp`` have been removed due to
their reliance on Python 2 and lack of maintenance.

What This Tool Does
-------------------

This tool helps you run molecular dynamics simulations in **LAMMPS** when you have your molecular system set up in **AMBER** format.

**Typical workflow:**

1. Start with a molecular structure (PDB file or SMILES string)
2. Use AMBER tools to create AMBER files (``.prmtop``, ``.mol2``, ``.frcmod``)
3. **Convert AMBER files to LAMMPS format** using this tool
4. Run your simulation in LAMMPS

**What gets converted:**

* **``.prmtop``** (topology file): Contains bonds, angles, atom types → LAMMPS data file
* **``.mol2``** (coordinates file): Contains atomic positions and charges → LAMMPS coordinates
* **``.frcmod``** (force field file): Contains interaction parameters → LAMMPS parameters

**What you get as output:**

* **LAMMPS data file** (e.g., ``data.lammps``): Contains atomic coordinates, box dimensions, and molecular topology
* **LAMMPS parameter file** (e.g., ``parm.lammps``): Contains force field parameters for bonds, angles, and nonbonded interactions



.. _amber2lammps: https://github.com/askforarun/AMBER2LAMMPS

Project and download
--------------------

AMBER2LAMMPS is developed and maintained outside the LAMMPS repository.
Download it from the upstream project and clone from GitHub:

.. code-block:: bash

   git clone https://github.com/askforarun/AMBER2LAMMPS.git
   cd AMBER2LAMMPS

Requirements
------------

**Platform Compatibility**

AMBER2LAMMPS works on all major platforms:

* **Linux**: Full support (Ubuntu, CentOS, Fedora, etc.)
* **macOS**: Full support (Intel and Apple Silicon)
* **Windows**: Support via WSL2 or Git Bash

**System Requirements**

* **Structure**: PDB file (or SMILES string, you can convert to PDB).
* **AmberTools utilities**: ``antechamber``, ``parmchk2``, ``tleap`` to
  build ``.prmtop``, ``.mol2``, and ``.frcmod``.
* **Python packages**: ``parmed`` and ``numpy``.
* **LAMMPS**: On your ``PATH`` (``which lmp`` or ``lmp -help``) and built
  with ``MOLECULE``, ``KSPACE``, and ``EXTRA-MOLECULE`` packages.
* **Optional**: Open Babel (``obabel``) if starting from SMILES.

Install dependencies
--------------------

AmberTools
^^^^^^^^^^

Install from https://ambermd.org/GetAmber.php#ambertools and activate
the environment:

.. code-block:: bash

   conda activate Ambertools23  # or your AmberTools environment

Python packages
^^^^^^^^^^^^^^^

Using conda (recommended):

.. code-block:: bash

   conda install -c conda-forge parmed numpy

Using pip:

.. code-block:: bash

   pip install parmed numpy

Open Babel (optional, SMILES → PDB)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using conda (recommended):

.. code-block:: bash

   conda install -c conda-forge openbabel

Using pip:

.. code-block:: bash

   pip install openbabel

Using system package manager:

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get install openbabel

   # macOS
   brew install open-babel

Command Reference
-----------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Argument
     - Required?
     - Description
   * - **Input Files**
     -
     -
   * - ``topology``
     - yes
     - AMBER topology file (``.prmtop``)
   * - ``mol2``
     - yes
     - MOL2 coordinate file (with charges)
   * - ``frcmod``
     - yes
     - AMBER force-field parameter file
   * - **Output Files**
     -
     -
   * - ``data_file``
     - yes
     - Output LAMMPS data filename
   * - ``param_file``
     - yes
     - Output LAMMPS parameter filename
   * - **Options**
     -
     -
   * - ``-b, --buffer``
     - optional
     - Vacuum padding (Å) for the simulation box. Default: ``3.8``.
   * - ``--verbose``
     - optional
     - Print step-by-step progress, counts, and box size. Default: ``False``.
   * - ``-h, --help``
     - optional
     - Show help message.

Conversion Process
------------------

The converter implements the following sequence (see
``amber_to_lammps.py``):

1. **Input validation**: ``validate_files`` confirms that ``topology``,
   ``mol2``, and ``frcmod`` exist and are readable. The CLI performs this validation automatically, while it's optional for the Python API.
2. **AMBER topology load**: ParmEd reads atoms, bonds, angles, and
   dihedrals from ``.prmtop``; counts are printed with ``--verbose``.
3. **Atom typing and masses**: ``MASS`` entries in ``frcmod`` map atom
   types to sequential IDs and masses; missing types warn and fall back
   to type 1.
4. **Coordinates and charges**: Coordinates and charges are read from the
   MOL2 ``ATOM`` block.
5. **Box creation (``--buffer``)**: A bounding box around the coordinates
   is expanded by the buffer on all sides.
6. **Charge normalization**: If the net charge magnitude exceeds
   ``1e-6``, charges are shifted uniformly to neutralize the system.
7. **Nonbonded parameters**: ``NONBON`` terms in ``frcmod`` become
   ``pair_coeff`` entries in the parameter file.
8. **Topology terms**: Bonds, angles, and dihedrals are exported via
   ParmEd to temporary files and written to the LAMMPS data/parameter
   files.
9. **Cleanup**: Temporary helper files (``bonds.txt``, ``angles.txt``,
   ``dihedrals.txt``) are removed.
10. **Verbose diagnostics**: With ``--verbose``, the script reports
    counts, box extents after buffering, and the normalized total charge.

Workflow Examples
-----------------

Prepare input files (ethanol example)
-------------------------------------

Convert SMILES to PDB (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you start from a SMILES string, generate a PDB with ``obabel``:

.. code-block:: bash

   obabel -:CCO -h -opdb -O ethanol.pdb --gen3d  # adds hydrogens, recommended
   obabel -:c1ccccc1 -h -opdb -O benzene.pdb --gen3d
   obabel -:"CC(=O)OC1=CC=CC=C1C(=O)O" -h -opdb -O aspirin.pdb --gen3d

Generate AMBER topology and parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assuming ``ethanol.pdb`` is your structure:

1. Generate a MOL2 file with charges:

   .. code-block:: bash

      antechamber -j 4 -at gaff2 -dr yes -fi pdb -fo mol2 \
         -i ethanol.pdb -o ethanol.mol2 -c bcc

   The ``-c bcc`` assigns AM1-BCC charges; include hydrogens to avoid missing
   types/charges.

2. Build the force-field modifications file:

   .. code-block:: bash

      parmchk2 -i ethanol.mol2 -o ethanol.frcmod -f mol2 -a Y

   The ``-a Y`` option is necessary to ensure all force-field parameters are included.

3. Create ``tleap.in``:

   .. code-block:: text

      source leaprc.gaff2
      SUS = loadmol2 ethanol.mol2
      check SUS
      loadamberparams ethanol.frcmod
      saveamberparm SUS ethanol.prmtop
      quit

4. Run tleap and inspect ``leap.log`` for errors:

   .. code-block:: bash

      tleap -f tleap.in

Outputs: ``ethanol.prmtop``, ``ethanol.mol2``, and
``ethanol.frcmod``.

Basic conversion workflow
--------------------------

LAMMPS input script
^^^^^^^^^^^^^^^^^^^^

Save the following LAMMPS commands in a file called ``example_lammps_input.lmp``:

.. code-block:: text

   # LAMMPS Input Script for Converted AMBER System

   units real
   dimension 3
   boundary p p p
   atom_style full

   read_data data.lammps

   pair_style      lj/cut/coul/long 9 9
   bond_style      harmonic
   angle_style     harmonic
   dihedral_style  fourier
   special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.83333333

   include parm.lammps

   thermo_style custom ebond eangle edihed eimp epair evdwl ecoul elong etail pe
   run 0

CLI usage with LAMMPS execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ensure amber_to_lammps.py is executable and the input files are accessible.

Run the converter:

.. code-block:: bash

   python3 amber_to_lammps.py data.lammps parm.lammps \
      ethanol.prmtop ethanol.mol2 ethanol.frcmod --verbose -b 4.5

Outputs: ``data.lammps`` (coordinates/topology) and ``parm.lammps``
(parameters). ``--buffer`` sets box padding in Å.

Run with LAMMPS:

.. code-block:: bash

   lmp < example_lammps_input.lmp

Additional CLI examples
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Custom buffer and verbose logging (adds 5 Å padding)
   python3 amber_to_lammps.py my_data.lammps my_params.lammps \
      ethanol.prmtop ethanol.mol2 ethanol.frcmod --verbose -b 5.0

   # Custom output names
   python3 amber_to_lammps.py system.data system.parm system.prmtop \
      system.mol2 system.frcmod

   # Minimal output without verbose logging
   python3 amber_to_lammps.py small.data small.parm \
      ethanol.prmtop ethanol.mol2 ethanol.frcmod -b 3.0

   # Using absolute paths with custom buffer
   python3 amber_to_lammps.py /home/user/lammps/output/data.lammps /home/user/lammps/output/param.lammps /home/user/amber/topology.prmtop /home/user/amber/coords.mol2 /home/user/amber/params.frcmod -b 4.5

Python API usage with LAMMPS execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from amber_to_lammps import amber2lammps, validate_files
   import subprocess

   validate_files('ethanol.prmtop', 'ethanol.mol2', 'ethanol.frcmod')

   amber2lammps(
       data_file='data.lammps',
       param_file='parm.lammps',
       topology='ethanol.prmtop',
       mol2='ethanol.mol2',
       frcmod='ethanol.frcmod',
       buffer=3.8,
       verbose=True,
   )

   cmd = "lmp < example_lammps_input.lmp"
   subprocess.call(cmd, shell=True)

Additional API examples
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Example 1: Basic conversion without validation
   from amber_to_lammps import amber2lammps

   amber2lammps(
       data_file='system.data',
       param_file='system.parm',
       topology='system.prmtop',
       mol2='system.mol2',
       frcmod='system.frcmod'
   )

   # Example 2: Custom buffer and verbose output
   amber2lammps(
       data_file='molecule.data',
       param_file='molecule.parm',
       topology='molecule.prmtop',
       mol2='molecule.mol2',
       frcmod='molecule.frcmod',
       buffer=5.0,
       verbose=True
   )

   # Example 3: Batch processing multiple molecules
   from amber_to_lammps import amber2lammps, validate_files

   molecules = ['ethanol', 'benzene', 'aspirin']

   for mol in molecules:
       validate_files(f'{mol}.prmtop', f'{mol}.mol2', f'{mol}.frcmod')
       amber2lammps(
           data_file=f'{mol}.data',
           param_file=f'{mol}.parm',
           topology=f'{mol}.prmtop',
           mol2=f'{mol}.mol2',
           frcmod=f'{mol}.frcmod',
           verbose=True
       )



Common issues and solutions
---------------------------

* **Missing MASS/atom types**: Check ``.frcmod`` and ensure -a Y option was used with ``parmchk2``.
* **Net charge not zero**: Charge normalization shifts charges; verify source charges and rerun if unintended.
* **Atoms too close or outside box**: Increase ``--buffer`` and rerun; inspect verbose box extents.
* **LAMMPS run errors about packages**: Rebuild LAMMPS with ``MOLECULE``, ``KSPACE``, ``EXTRA-MOLECULE``.

Validation
----------

AMBER2LAMMPS has been validated against InterMol output. See the project
page for details: https://github.com/askforarun/AMBER2LAMMPS

Getting Help
------------

* **Submit Issues:** https://github.com/askforarun/AMBER2LAMMPS/issues
* **Feature Requests:** Use GitHub Issues or Discussions
* **Questions:** Use GitHub Discussions or Issues

Citation
--------

If you use this Python tool in your research, please cite it as:

**DOI:** 10.5281/zenodo.18114886

License
-------

AMBER2LAMMPS is distributed under the MIT license.

.. raw:: latex

   \clearpage
