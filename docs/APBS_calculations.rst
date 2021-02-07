.. _APBS:

Running APBS calculations
=========================

APBS input Generation
---------------------

After :ref:`System Setup <system_setup>`, APBS calculations are performed to
obtain the relevant information. The input of the APBS calculations can be
generated either automatically or manually.

.. toctree::
   :maxdepth: 1
   :caption: Generation of APBS Input

   automatic_APBS
   manual_APBS

.. _run_APBS:

Run APBS calculations and load the result
-----------------------------------------

After the generation of the APBS input file (:code:`apbs_in='apbs.in'`),
the APBS calculation could be perfomed with
:func:`rocklinc.RocklinCorrection.run_APBS` and the result is loaded to the
object via :func:`rocklinc.RocklinCorrection.read_APBS`.

For automatic APBS calculation:

.. code-block:: python

    import rocklinc
    rc = rocklinc.RocklinCorrection(box, lig_netq, protein_netq, temp,
                                    rocklinc.waters.TIP3P)
    rc.make_APBS_input(u, 'resname LIG',)
    rc.run_APBS()
    rc.read_APBS()

For manual APBS calculation:

.. code-block:: python

    import rocklinc
    rc = rocklinc.RocklinCorrection(box, lig_netq, protein_netq, temp,
                                    rocklinc.waters.TIP3P)
    rc.set_APBS_input(800, box, lig_netq, protein_netq,)
    rc.run_APBS()
    rc.read_APBS()

.. autofunction:: rocklinc.RocklinCorrection.run_APBS
.. autofunction:: rocklinc.RocklinCorrection.read_APBS

After loading the results from APBS calculation, the result is
:ref:`summarised <summary>` and can be written into a text file.