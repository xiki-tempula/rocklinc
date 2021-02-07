.. _summary:

Compiling Results
=================

After the APBS calculations, the result could be return via the
:func:`rocklinc.RocklinCorrection.compute` and the detailed break down could
be written to a text file with :func:`rocklinc.RocklinCorrection.write`.

For automatic APBS calculation:

.. code-block:: python

    import rocklinc
    rc = rocklinc.RocklinCorrection(box, lig_netq, protein_netq, temp,
                                    rocklinc.waters.TIP3P)
    rc.make_APBS_input(u, 'resname LIG',)
    rc.run_APBS()
    rc.read_APBS()
    result = rc.compute()
    rc.write('correction.txt')

`result` is the correction energy and `'correction.txt'` has the format::

    The total correction energy is: 2.62 kJ/mol or 0.63 kCal/mol
    = ΔΔG_ANA(L): 2.62 kJ/mol + ΔΔG_DSC(L): 0.00 kJ/mol
    ΔΔG_ANA(L) = ΔΔG_NET(L) + ΔΔG_USV(L) + ΔΔG_RIP(L) + ΔΔG_EMP(L)
    ΔΔG_NET(L) = 98.55 kJ/mol
    ΔΔG_NET(L) + ΔΔG_USV(L) = 1.02 kJ/mol
    ΔΔG_RIP(L) = 1.61 kJ/mol
    ΔΔG_EMP(L) = -0.00 kJ/mol

.. autofunction:: rocklinc.RocklinCorrection.compute
.. autofunction:: rocklinc.RocklinCorrection.write

