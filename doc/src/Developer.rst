Information for Developers
==========================

This section describes the internal structure and basic algorithms
of the LAMMPS code. This is a work in progress and additional
information will be added incrementally depending on availability
of time and requests from the LAMMPS user community.

A discussion of software engineering methods applied to LAMMPS over time
and a general outline of the design and maintenance approach by the
LAMMPS developers can be found in the paper `LAMMPS: A Case Study For
Applying Modern Software Engineering to an Established Research Software
Package, <https://doi.org/10.5281/zenodo.17117558>`_ in the USRSE'25
conference proceedings.

.. admonition:: On using AI coding agents
   :class: Hint

   Recently AI coding agents have seen massive improvements and thus
   have become quite popular and viable in assisting with the
   development of new code or with adding features to existing code.
   The LAMMPS developers have made some experiments and tried to
   assess the how effective and reliable the results are.  Below are
   some of our observations.

   A significant part of LAMMPS was written by taking some existing
   classes that perform similar calculations or operations and then
   replacing the functionality.  That kind of task is well suited for a
   coding agent, but there are often little details that require some
   manual adjustments and thus careful checking of the resulting code is
   recommended.  Similarly, porting styles to the KOKKOS or OPENMP
   package can be delegated in many cases, provided suitable template
   examples exist and the coding agent is asked to look at those.
   Asking a coding agent to create new functionality for which there
   is no similar code in either the LAMMPS sources or some other
   repository has a much lower chance of success.

   Another task that can often be passed successfully to a coding agent
   is the creation of suitable documentation.  However, the results have
   a tendency to be overly verbose and sometimes redundant, so some
   manual editing to tighten it is recommended.

   Finally, coding agents are quite good at reviewing existing code.
   The LAMMPS developers use GitHub Copilot regularly to review
   submitted pull requests and over time the rate of incorrect advice
   has become rather low.  When using an AI coding agent, it is
   generally recommended to use a different model for the code review
   than what was used to write the code.

   To assist coding agents in their tasks, and particularly to provide
   instructions that GitHub Copilot will automatically read, we have
   created a file ``.github/copilot-instructions.md`` that can be used
   to inform coding agents about important steps when creating or
   reviewing code for LAMMPS.

   When contributing code to LAMMPS that was at least partially created
   by a coding agent, please provide information in the pull request
   description to the extent of AI usage and which model or models were
   used.

.. toctree::
   :maxdepth: 1

   Developer_org
   Developer_code_design
   Developer_parallel
   Developer_atom
   Developer_comm_ops
   Developer_flow
   Developer_write
   Developer_notes
   Developer_updating
   Developer_plugins
   Developer_unittest
   Classes
   Developer_platform
   Developer_utils
   Developer_internal
   Developer_grid
