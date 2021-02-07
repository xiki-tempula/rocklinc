Automatic Generation of APBS Input Files
========================================

The APBS input can be generated automatically with :class:`MDAnalysis.Universe`
object, which should contain the attribute `positions`, `charges`, `radii` and
`dimensions`. Within the :class:`MDAnalysis.Universe`, different molecules
should be able to be accessed through the `residues` attribute.

Tips for Universe generation
----------------------------

One way of generating such a universe via Gromacs is to obtain `charges` from
the `tpr` file, the `positions` from `gro` file and `radii` from `pqr` files.

The `pqr` files could be generated using `editconf` functionality.

.. code-block:: bash

    gmx editconf -f topol.tpr -mead topol.pqr

The :class:`MDAnalysis.Universe` could then be assembled:

.. code-block:: python

    import MDAnalysis as mda
    u = mda.Universe('topol.tpr', 'topol.gro')
    pqr = mda.Universe('topol.pqr')
    u.add_TopologyAttr('radii')
    u.atoms.radii = pqr.atoms.radii

.. note::
    Only the pqr files generated with certain versions of Gromacs could be
    loaded by MDAnalysis and `Gromacs 5.1.5` is one of them.

APBS Input Generation
---------------------
The APBS input could then be generated with
:func:`rocklinc.RocklinCorrection.make_APBS_input`.

.. code-block:: python

    import rocklinc
    rc = rocklinc.RocklinCorrection(box, lig_netq, protein_netq, temp,
                                    rocklinc.waters.TIP3P)
    rc.make_APBS_input(u, 'resname LIG',)

The :func:`rocklinc.RocklinCorrection.make_APBS_input` also sets the variable
:attr:`rocklinc.RocklinCorrection.NS`, which is the number of solvent molecules.
If there isn't any none solvent and ligand molecule, the integrated potential of
the protein :attr:`rocklinc.RocklinCorrection.IP` will be set to 0 as well.

.. autofunction:: rocklinc.RocklinCorrection.make_APBS_input

The next step is to :ref:`run the APBS calculations <run_APBS>`.