.. index:: fix mbx

    fix mbx command
    ===============

    Syntax
    """"""

    .. code-block:: LAMMPS
        fix ID group-ID mbx 

* ID, group-ID are documented in :doc:`fix <fix>` command








Restrictions
""""""""""""

This fix is part of the MBX package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info. This fix also relies on the
presence of pair_mbx command

There can only be one fix mbx command active at a time.




Related commands
""""""""""""""""

:doc:`pair mbx <pair_mbx>`

Default
"""""""


-----------