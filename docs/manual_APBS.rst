Manual Generation of APBS Input Files
=====================================

The APBS input could be set manually as well. Three `pqr` files (
`prot_only.pqr`, `lig_in_prot.pqr`, `lig_only.pqr`) has to be generated
according to the original definition from . For a system that has 800 water
molecules, the files could be feed to the object via
:func:`rocklinc.RocklinCorrection.set_APBS_input`

.. code-block:: python

    import rocklinc
    rc = rocklinc.RocklinCorrection(box, lig_netq, protein_netq, temp,
                                    rocklinc.waters.TIP3P)
    rc.set_APBS_input(800, box, lig_netq, protein_netq,)

.. autofunction:: rocklinc.RocklinCorrection.set_APBS_input

The next step is to :ref:`run the APBS calculations <run_APBS>`.