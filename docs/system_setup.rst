.. _system_setup:

System Setup
===============

Prepare Input
-------------
Aside from the APBS calculation, Rocklin correction required several input

* box: The size of the simulation box is the average box dimensions across the
  simulations. The input is a list with length 3 (e.g. `[10, 10, 10]`). The
  unit should be Å.
* lig_netq: The charge of the ligand.
* protein_netq: The charge of the system without the ligand.
* temp: The temperature of the system in kelvin (e.g. 310K)
* water: The water model used in the simulation. Currently only tip3p and tip4p
  water are supported. The input should be :mod:`rocklinc.waters` objects.

Run command
-----------
To use `rocklinc`, one need to create a :class:`rocklinc.RocklinCorrection`
object to perform the computation.

.. code-block:: python

    import rocklinc
    rc = rocklinc.RocklinCorrection(box, lig_netq, protein_netq, temp,
                                    rocklinc.waters.TIP3P)

This is followed up by :ref:`APBS calculations <APBS>`.

Optional Water models
---------------------

.. automodule:: rocklinc.waters
   :members: TIP3P, TIP4P

