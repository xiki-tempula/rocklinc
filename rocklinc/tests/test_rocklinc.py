"""
Unit and regression test for the rocklinc package.
"""

# Import package, test suite, and other packages as needed
import rocklinc
import pytest
import MDAnalysis as mda
import sys
from io import StringIO
import quantities as pq
from pkg_resources import resource_filename
import numpy as np
import os

def test_manual_cl():
    box = [2, 2, 2] * pq.nm
    lig_netq = -1 * pq.e
    protein_netq = 0 * pq.e
    temp = 310 * pq.kelvin
    new = rocklinc.RocklinCorrection(box, lig_netq, protein_netq, temp)

    new.set_APBS_input(0)
    # check apbs.in
    with open('apbs.in', 'r') as f:
        current_apbs_in = f.read()
    with open(resource_filename(__name__, 'test_CL/apbs.in'), 'r') as f:
        model_apbs_in = f.read()
    assert current_apbs_in == model_apbs_in



def test_automatic_cl():
    gro='''Single cl molecule                                                     
1                                                                               
    1  CL    CL    1   1.000   1.000   1.000                   
2 2 2  '''
    u = mda.Universe(StringIO(gro), format='gro')
    u.add_TopologyAttr('charges')
    u.select_atoms('resname CL').charges = -1
    u.add_TopologyAttr('radii')
    u.select_atoms('resname CL').radii = 2.20

    box = [2, 2, 2] * pq.nm
    lig_netq = -1 * pq.e
    protein_netq = 0 * pq.e
    temp = 310 * pq.kelvin
    lig_selection = 'resname CL'

    new = rocklinc.RocklinCorrection(box, lig_netq, protein_netq, temp)
    new.make_APBS_input(u.select_atoms('resname SOL CL NA'), lig_selection,)

    # check apbs.in
    with open('apbs.in', 'r') as f:
        current_apbs_in = f.read()
    with open(resource_filename(__name__, 'test_CL/apbs.in'), 'r') as f:
        model_apbs_in = f.read()
    assert current_apbs_in == model_apbs_in

    # check lig_in_prot.pqr
    with open('lig_in_prot.pqr', 'r') as f:
        current_lig_in_prot = f.read().strip().split('\n')[-2:]
    with open(resource_filename(__name__, 'test_CL/lig_in_prot.pqr'), 'r') as f:
        model_lig_in_prot = f.read().strip().split('\n')[-2:]
    assert current_lig_in_prot == model_lig_in_prot

    # check lig_only.pqr
    with open('lig_only.pqr', 'r') as f:
        current_lig_only = f.read().strip().split('\n')[-2:]
    with open(resource_filename(__name__, 'test_CL/lig_only.pqr'), 'r') as f:
        model_lig_only = f.read().strip().split('\n')[-2:]
    assert current_lig_only == model_lig_only

    # check prot_only.pqr
    with open('prot_only.pqr', 'r') as f:
        current_prot_only = f.read().strip().split('\n')[-2:]
    with open(resource_filename(__name__, 'test_CL/prot_only.pqr'), 'r') as f:
        model_prot_only = f.read().strip().split('\n')[-2:]
    assert current_prot_only == model_prot_only

    # new.run_APBS()

    new.read_APBS(ligand_RIP_het=resource_filename(__name__, 'test_CL/ligand_RIP_het.dx'),
                  protein_RIP_het=resource_filename(__name__, 'test_CL/protein_RIP_het.dx'),
                  ligand_RIP_hom=resource_filename(__name__, 'test_CL/ligand_RIP_hom.dx'))
    result = new.compute()
    os.remove('prot_only.pqr')
    os.remove('lig_only.pqr')
    os.remove('lig_in_prot.pqr')
    os.remove('apbs.in')
    np.testing.assert_almost_equal(result.base, 626.979, decimal=3)

    new.write('correction.txt')
    # check correction.txt
    with open('correction.txt', 'r') as f:
        current_correction = f.read()
    with open(resource_filename(__name__, 'test_CL/correction.txt'), 'r') as f:
        model_correction = f.read()
    assert current_correction == model_correction
    os.remove('correction.txt')