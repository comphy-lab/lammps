.. _howto_amber2lammps:
.. _howto_amber_to_lammps:

How-to: AMBER2LAMMPS
====================

AMBER2LAMMPS is a modern Python utility for converting AMBER molecular dynamics files to LAMMPS data format, replacing the legacy scripts in the ``amber2lmp`` folder shipped with LAMMPS with enhanced features and improved usability.

**Original Project:** https://github.com/askforarun/AMBER2LAMMPS

Overview
--------

AMBER2LAMMPS converts AMBER topology files (.prmtop), MOL2 coordinate files, and force field parameter files (.frcmod) into LAMMPS data and parameter formats. This tool provides enhanced features including charge normalization, separate output files, and improved error handling compared to the legacy scripts in the ``amber2lmp`` folder shipped with LAMMPS.

Key improvements over the legacy tool:

* **Charge normalization** - Automatically normalizes atomic charges to ensure zero net charge
* **Separate data and parameter files** - Generates distinct LAMMPS data file and parameter file for better organization
* **Professional CLI** - Comprehensive command-line interface with proper error handling
* **Verbose output** - Detailed debugging and monitoring information
* **Cross-platform** - Works on Linux, macOS, and Windows

Installation
------------

Prerequisites
"""""""""""""

* Python 3.6 or higher
* pip package manager

Required Dependencies
"""""""""""""""""""""

Install the required packages using pip (recommended):

.. code-block:: bash

   pip install numpy parmed

Or using conda:

.. code-block:: bash

   conda install numpy
   pip install parmed

Obtaining AMBER2LAMMPS
""""""""""""""""""""""

AMBER2LAMMPS is developed and maintained outside the LAMMPS repository.
Download it from the upstream project:

.. code-block:: bash

   git clone https://github.com/askforarun/AMBER2LAMMPS.git
   cd AMBER2LAMMPS

AmberTools (Optional)
"""""""""""""""""""""

If you have AMBERTools installed, you can activate the environment:

.. code-block:: bash

   conda activate Ambertools23  # or your AMBERTools version

Platform Compatibility
----------------------

AMBER2LAMMPS has been validated and tested on:

* **Linux** (Ubuntu, CentOS, Red Hat, Debian) - Fully tested and validated
* **macOS** (Intel and Apple Silicon) - Fully tested
* **Windows** (Windows 10/11 with WSL2 and native Python) - Tested with WSL2 and native Python

WSL2 is recommended for Windows users for best compatibility.

Usage
-----

Command Line Interface
""""""""""""""""""""""

Basic Usage
^^^^^^^^^^^

The following examples assume you are running from an AMBER2LAMMPS checkout
or otherwise have access to its ``amber_to_lammps.py`` script.

.. code-block:: bash

   python3 amber_to_lammps.py <data_file> <param_file> <topology> <mol2> <frcmod>

Arguments
^^^^^^^^^

**Positional Arguments:**

* ``data_file``: Output LAMMPS data file name
* ``param_file``: Output LAMMPS parameter file name  
* ``topology``: AMBER topology file (.prmtop)
* ``mol2``: MOL2 coordinate file
* ``frcmod``: Force field parameter file (.frcmod)

**Optional Arguments:**

* ``-b, --buffer``: Buffer size around molecule (default: 3.8)
* ``--verbose``: Enable verbose output for debugging

Example
^^^^^^^

.. code-block:: bash

   python3 amber_to_lammps.py system.data system.param system.prmtop system.mol2 system.frcmod -b 4.0 --verbose

Python API
""""""""""

Basic Usage
^^^^^^^^^^^

.. code-block:: python

   from amber_to_lammps import amber2lammps

   amber2lammps(
       data_file="system.data",
       param_file="system.param", 
       topology="system.prmtop",
       mol2="system.mol2",
       frcmod="system.frcmod",
       buffer=3.8,
       verbose=True
   )

Advanced Usage
^^^^^^^^^^^^^^

.. code-block:: python

   from amber_to_lammps import amber2lammps
   import os

   # Input files
   topology = "system.prmtop"
   mol2 = "system.mol2"
   frcmod = "system.frcmod"

   # Output files
   data_file = "system.data"
   param_file = "system.param"

   # Convert with custom buffer and verbose output
   amber2lammps(
       data_file=data_file,
       param_file=param_file,
       topology=topology,
       mol2=mol2,
       frcmod=frcmod,
       buffer=4.0,
       verbose=True
   )

   # Check if conversion was successful
   if os.path.exists(data_file) and os.path.exists(param_file):
       print("Conversion completed successfully!")
       print(f"Output files: {data_file}, {param_file}")
   else:
       print("Conversion failed!")

Generating AMBER Input Files
---------------------------

If you need to generate AMBER input files from PDB, use this workflow:

.. code-block:: python

   import subprocess

   # Step 1: Generate MOL2 file with charges
   cmd1 = "antechamber -j 4 -at gaff2 -dr no -fi pdb -fo mol2 -i molecule.pdb -o molecule.mol2 -c bcc"
   subprocess.run(cmd1, shell=True)

   # Step 2: Generate force field parameters
   cmd2 = "parmchk2 -i molecule.mol2 -o molecule.frcmod -f mol2 -a Y"
   subprocess.run(cmd2, shell=True)

   # Step 3: Create tleap input file
   with open("tleap.in", "w") as f:
       f.write("source leaprc.gaff2\\n")
       f.write("MOL = loadmol2 molecule.mol2\\n") 
       f.write("check MOL\\n")
       f.write("loadamberparams molecule.frcmod\\n")
       f.write("saveamberparm MOL molecule.prmtop molecule.crd\\n")
       f.write("quit")

   # Step 4: Run tleap to generate AMBER files
   cmd3 = "tleap -f tleap.in"
   subprocess.run(cmd3, shell=True)

**Prerequisites for AMBER Workflow:**

* AMBERTools must be installed and in PATH
* Input PDB file: ``molecule.pdb``
* Output files: ``molecule.prmtop``, ``molecule.crd``, ``molecule.mol2``, ``molecule.frcmod``

Using in LAMMPS
---------------

After conversion, you can use the generated files in your LAMMPS simulation:

.. code-block:: lammps

   # LAMMPS Input Script for Converted AMBER System

   # --------------------- Initialization ---------------------
   units real
   dimension 3
   boundary p p p
   atom_style full

   # --------------------- System Setup ---------------------
   # Read the molecular structure
   read_data data.lammps

   # --------------------- Force Field Settings ---------------------
   pair_style      lj/cut/coul/long 9 9
   kspace_style    pppm 1.0e-8
   pair_modify     tail yes
   bond_style      harmonic
   angle_style      harmonic
   dihedral_style    fourier
   special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.83333333

   # Include the force field parameters generated by amber_to_lammps
   include parm.lammps

   # --------------------- Output Settings ---------------------
   thermo_style custom ebond eangle edihed eimp epair evdwl ecoul elong etail pe

   # --------------------- Energy Test ---------------------
   # Run 0 steps to check energy and system setup
   run 0

Troubleshooting
---------------

Common Issues
"""""""""""""

**Error: "Cannot find topology file"**

* Ensure the .prmtop file exists and is readable
* Check file permissions
* Verify the file path is correct

**Error: "Cannot find MOL2 file"**

* Ensure the .mol2 file exists and is readable
* Check that the MOL2 file was generated correctly by AMBERTools
* Verify the file path is correct

**Error: "Cannot find frcmod file"**

* Ensure the .frcmod file exists and is readable
* Check that parmchk2 generated the file correctly
* Verify the file path is correct

**Charges not normalized properly**

* Check the input topology file for correct charge information
* Verify that the MOL2 file contains charge information
* Use the ``--verbose`` flag to debug charge normalization

Validation
""""""""""

AMBER2LAMMPS has been validated against InterMol output for accuracy:

* **Charge normalization** - Ensures zero net charge
* **Bond parameters** - Correctly extracted from AMBER force field
* **Angle parameters** - Properly converted to LAMMPS format
* **Dihedral parameters** - Accurately translated from AMBER

Getting Help
""""""""""""

For additional help, bug reports, feature requests, or questions about AMBER2LAMMPS:

* **Submit Issues:** https://github.com/askforarun/AMBER2LAMMPS/issues
* **Feature Requests:** Use GitHub Issues or Discussions
* **Questions:** Use GitHub Discussions or Issues

Citation
--------

If you use this software in your research, please cite it as:

**DOI:** 10.5281/zenodo.18114886

License
-------

AMBER2LAMMPS is MIT licensed (external).
