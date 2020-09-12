"""
rocklin_correction.py
A python module for performing Rocklin Correction.

Handles the primary functions
"""

from pkg_resources import resource_filename
from subprocess import call

import quantities as pq
import numpy as np
from gridData import Grid
import MDAnalysis as mda

from . import constants
from .waters import TIP3P


class RocklinCorrection():
    def __init__(self, box, lig_netq, protein_netq, temp=None, water=None):
        try:
            box.rescale(pq.angstrom)
            self.box = box
        except AttributeError:
            self.box = box * pq.angstrom
        self.vol = self.box[0] * self.box[1] * self.box[2]
        try:
            lig_netq.rescale(pq.e)
            self.lig_netq = lig_netq
        except AttributeError:
            self.lig_netq = lig_netq * pq.e

        try:
            protein_netq.rescale(pq.e)
            self.protein_netq = protein_netq
        except AttributeError:
            self.protein_netq = protein_netq * pq.e

        if temp is None:
            self.temp = 298.15 * pq.Kelvin
        else:
            try:
                temp.rescale(pq.Kelvin)
                self.temp = temp
            except AttributeError:
                self.temp = temp * pq.Kelvin
        if water is None:
            self.water = TIP3P
        else:
            self.water = water

    def make_APBS_input(self, tpr, cord, pqr, ligand_selection,
                        solvent_selection='resname SOL',
                        out_prot_only='prot_only.pqr',
                        out_lig_in_prot='lig_in_prot.pqr',
                        out_lig_only='lig_only.pqr'):
        self.out_prot_only = out_prot_only
        self.out_lig_in_prot = out_lig_in_prot
        self.out_lig_only = out_lig_only


        box = cord.dimensions[:3]
        self.apbs_box = cord.dimensions[:3] * pq.angstrom
        self.apbs_vol = np.prod(self.apbs_box)
        pqr.atoms.positions = cord.atoms.positions
        pqr.dimensions = cord.dimensions

        # Charge only the ligand
        pqr.atoms.charges = tpr.atoms.charges
        pqr.select_atoms('not {}'.format(ligand_selection)).charges = 0
        pqr.select_atoms('not {}'.format(solvent_selection)).write(out_lig_in_prot)
        # Charge only the Rest of the system
        pqr.atoms.charges = tpr.atoms.charges
        pqr.select_atoms('{}'.format(ligand_selection)).charges = 0
        pqr.select_atoms('not {}'.format(solvent_selection)).write(out_prot_only)

        # Ligand for centering
        pqr.select_atoms('{}'.format(ligand_selection)).write(out_lig_only)
        self.NS = len(tpr.select_atoms(solvent_selection).residues)

    def run_APBS(self, apbs_exe='/opt/local/bin/apbs', apbs_in='apbs.in', apbs_out='apbs.out', box=None):
        with open(resource_filename(__name__, 'data/apbs.in'), 'r') as f:
            txt = f.read()
        if box is None:
            box = self.apbs_box
        with open(apbs_in, 'w') as f:
            f.write(txt.format(prot_only=self.out_prot_only,
                               lig_in_prot=self.out_lig_in_prot,
                               lig_only=self.out_lig_only,
                               x=box[0].base, y=box[1].base, z=box[2].base,
                               e=self.water.epsilon_S,
                               t=self.temp.base,))
        call([apbs_exe, apbs_in])

    def dx2IP(self, dx):
        g = Grid(dx)
        V = ((np.prod(g.delta*pq.angstrom)) * np.prod(g.grid.shape))
        return np.average(g.grid) * (1/pq.e) * (constants.kB * self.temp) * V

    def read_APBS(self, ligand_RIP_het='ligand_RIP_het.dx',
                  protein_RIP_het='protein_RIP_het.dx',
                  ligand_RIP_hom='ligand_RIP_hom.dx'):
        # Ligand Het
        IL_Bx = self.dx2IP(ligand_RIP_het)
        IL_BQx = (-constants.xi_CB * constants.coulomb_factor / self.water.epsilon_S) * self.lig_netq * (self.apbs_vol ** (2.0 / 3.0))
        self.IL = IL_Bx - IL_BQx
        # Protein Het
        IP_Bx = self.dx2IP(protein_RIP_het)
        IP_BQx = (-constants.xi_CB * constants.coulomb_factor / self.water.epsilon_S) * self.protein_netq * (self.apbs_vol ** (2.0 / 3.0))
        self.IP = IP_Bx - IP_BQx
        # Ligand Het
        IL_hom_Bx = self.dx2IP(ligand_RIP_hom)
        IL_hom_BQx = (-constants.xi_CB * constants.coulomb_factor / 1) * self.lig_netq * (self.apbs_vol ** (2.0 / 3.0))
        IL_hom = IL_hom_Bx - IL_hom_BQx
        self.IL_SLV = self.IL - IL_hom

    def compute(self):
        delta_DSC = - (self.water.gamma_s * self.lig_netq) / (6 * constants.epsilon_0) * self.NS / self.vol
        delta_NET = - constants.xi_LS / (8 * np.pi * constants.epsilon_0) * (
                    (self.protein_netq + self.lig_netq) ** 2 - self.protein_netq ** 2) / (self.vol ** (1/3))
        delta_NET_delta_USV = delta_NET / self.water.epsilon_S
        delta_RIP = ((self.IP + self.IL) * (self.protein_netq + self.lig_netq) - self.IP * self.protein_netq) / self.apbs_vol
        RL = ((1 / (8 * np.pi * constants.epsilon_0) * (4 * np.pi / 3) * (
                    1 - 1 / self.water.epsilon_S) * self.lig_netq) ** -1 * self.IL_SLV) ** 0.5
        delta_EMP = - 1 / (8 * np.pi * constants.epsilon_0) * (16 * np.pi ** 2 / 45) * (
                    1 - 1 / self.water.epsilon_S) * ((self.protein_netq + self.lig_netq) ** 2 - self.protein_netq ** 2) * \
                    RL ** 5 / self.vol ** 2
        delta_ANA = delta_NET_delta_USV + delta_RIP + delta_EMP
        delta = delta_ANA + delta_DSC

        self.output = []
        self.output.append(
            'The total correction energy is: {:.2f} kJ/mol or {:.2f} kCal/mol'.format(
                delta.rescale(pq.J / pq.mol).item() / 1000,
                delta.rescale(pq.cal / pq.mol).item() / 1000))
        self.output.append(
            '= ΔΔG_ANA(L): {:.2f} kJ/mol + ΔΔG_DSC(L): {:.2f} kJ/mol'.format(
                delta_ANA.rescale(pq.J / pq.mol).item() / 1000,
                delta_DSC.rescale(pq.J / pq.mol).item() / 1000))
        self.output.append(
            'ΔΔG_ANA(L) = ΔΔG_NET(L) + ΔΔG_USV(L) + ΔΔG_RIP(L) + ΔΔG_EMP(L)')
        self.output.append('ΔΔG_NET(L) = {:.2f} kJ/mol'.format(
            delta_NET.rescale(pq.J / pq.mol).item() / 1000))
        self.output.append('ΔΔG_NET(L) + ΔΔG_USV(L) = {:.2f} kJ/mol'.format(
            delta_NET_delta_USV.rescale(pq.J / pq.mol).item() / 1000))
        self.output.append('ΔΔG_RIP(L) = {:.2f} kJ/mol'.format(
            delta_RIP.rescale(pq.J / pq.mol).item() / 1000))
        self.output.append('ΔΔG_EMP(L) = {:.2f} kJ/mol'.format(
            delta_EMP.rescale(pq.J / pq.mol).item() / 1000))
        return delta.rescale(pq.calorie / pq.mol)

    def write(self, outfile):
        with open(outfile, 'w') as f:
            f.write('\n'.join(self.output))


if __name__ == "__main__":
    tpr = mda.Universe('prod.tpr')
    pqr = mda.Universe('prod.pqr')
    cord = mda.Universe('prod.gro')
    lig_netq = -1 * pq.e
    protein_netq = 0 * pq.e
    temp = 310 * pq.kelvin
    lig_selection = 'resname CL'
    new = RocklinCorrection(cord.dimensions[:3], lig_netq, protein_netq,
                            temp)
    new.make_APBS_input(tpr, cord, pqr, lig_selection)
    new.run_APBS()
    new.read_APBS()
    result = new.compute()
    new.write('correction.txt')