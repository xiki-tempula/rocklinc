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
        '''
        Parameters
        ----------
        box : array
            The unitcell dimensions of the system ``[lx, ly, lz]``.
        lig_netq: float
            The unit charge of the ligand.
        protein_netq : float
            The unit charge of the protein.
        temp : float
            The temperature of the system in K.
        water : water
            The water model being used.
        '''
        self.box = pq.Quantity(box, pq.angstrom)
        self.vol = self.box[0] * self.box[1] * self.box[2]
        self.lig_netq = pq.Quantity(lig_netq, pq.e)
        self.protein_netq = pq.Quantity(protein_netq, pq.e)

        if temp is None:
            self.temp = 298.15 * pq.Kelvin
        else:
            self.temp = pq.Quantity(temp, pq.Kelvin)
        if water is None:
            self.water = TIP3P
        else:
            self.water = water

    def set_APBS_input(self, NS, box=None, qL=None, qP=None,
                       out_prot_only='prot_only.pqr',
                       out_lig_in_prot='lig_in_prot.pqr',
                       out_lig_only='lig_only.pqr',
                       apbs_in='apbs.in'):
        ''' Manually set the input file for the APBS calculations.
        Parameters
        ----------
        NS : int
            The number of solvent molecule in the system.
        box : array, optional
            The unitcell dimensions of the system ``[lx, ly, lz]`` for APBS
            calculations.
        qL: float, optional
            The unit charge of the ligand.
        qP : float, optional
            The unit charge of the protein.
        out_prot_only : str, optional
            The name of the pqr file where the ligand has no partial charge.
            (``prot_only.pqr``)
        out_lig_in_prot : str, optional
            The name of the pqr file where the protein has no partial charge.
            (``lig_in_prot.pqr``)
        out_lig_only : str, optional
            The name of the pqr file of the ligand.
            (``lig_only.pqr``)
        apbs_in: str, optional
            The input file to the APBS program. (``apbs.in``)
        '''
        self.NS = NS

        if box is None:
            self.apbs_box = self.box
        else:
            self.apbs_box = pq.Quantity(box, pq.angstrom)

        if qL is None:
            self.apbs_qL = self.lig_netq
        else:
            self.apbs_qL = pq.Quantity(qL, pq.e)

        if qP is None:
            self.apbs_qP = self.protein_netq
        else:
            self.apbs_qP = pq.Quantity(qP, pq.e)

        self.apbs_vol = np.prod(self.apbs_box)

        self.out_prot_only = out_prot_only
        self.out_lig_in_prot = out_lig_in_prot
        self.out_lig_only = out_lig_only
        self.IP = None

        self._write_APBS_input(apbs_in)


    def make_APBS_input(self, universe, ligand_selection,
                        solvent_selection='resname SOL',
                        out_prot_only='prot_only.pqr',
                        out_lig_in_prot='lig_in_prot.pqr',
                        out_lig_only='lig_only.pqr',
                        apbs_in='apbs.in'):
        ''' Automatically setup the input file for the APBS calculations.
        Parameters
        ----------
        universe : MDAnalysis.Universe
            The Universe object of the system, where the coordinate,
            partial charge, radii and system dimension will be obtained.
        ligand_selection: str
            The selection string for the ligand.
        solvent_selection : str, optional
            The selection string for the solvent. (``resname SOL``)
        out_prot_only : str, optional
            The name of the pqr file where the ligand has no partial charge.
            (``prot_only.pqr``)
        out_lig_in_prot : str, optional
            The name of the pqr file where the protein has no partial charge.
            (``lig_in_prot.pqr``)
        out_lig_only : str, optional
            The name of the pqr file of the ligand.
            (``lig_only.pqr``)
        apbs_in: str, optional
            The input file to the APBS program. (``apbs.in``)

        Attributes
        ----------
        IP : float
            The integrated potential of protein will be set to 0, if no
            protein is found.
        NS : int
            The number of solvent molecules in the system.
        '''
        self.out_prot_only = out_prot_only
        self.out_lig_in_prot = out_lig_in_prot
        self.out_lig_only = out_lig_only


        box = universe.dimensions[:3]
        self.apbs_box = box * pq.angstrom
        charges = universe.atoms.charges
        

        # Charge only the ligand
        universe.select_atoms('not {}'.format(ligand_selection)).charges = 0
        self.apbs_qL = np.sum(universe.atoms.charges) * pq.e
        universe.select_atoms('not {}'.format(solvent_selection)).write(out_lig_in_prot)
        # Charge only the Rest of the system
        universe.atoms.charges = charges
        universe.select_atoms('{}'.format(ligand_selection)).charges = 0
        self.apbs_qP = np.sum(universe.atoms.charges) * pq.e
        universe.select_atoms('not {}'.format(solvent_selection)).write(out_prot_only)

        # Check if there is anything other than ligand and solvent
        if len(universe.select_atoms('not (({}) or ({}))'.format(ligand_selection, solvent_selection))) > 0:
            self.IP = None
        else:
            self.IP = 0

        # Ligand for centering
        universe.select_atoms('{}'.format(ligand_selection)).write(out_lig_only)
        self.NS = len(universe.select_atoms(solvent_selection).residues)

        self._write_APBS_input(apbs_in)

    def _write_APBS_input(self, apbs_in):
        with open(resource_filename(__name__, 'data/apbs.in'), 'r') as f:
            txt = f.read()
        box = self.apbs_box
        with open(apbs_in, 'w') as f:
            f.write(txt.format(prot_only=self.out_prot_only,
                               lig_in_prot=self.out_lig_in_prot,
                               lig_only=self.out_lig_only,
                               x=box[0].magnitude, y=box[1].magnitude, z=box[2].magnitude,
                               e=self.water.epsilon_S,
                               t=self.temp.magnitude,))

    def run_APBS(self, apbs_exe='/opt/local/bin/apbs', apbs_in='apbs.in'):
        ''' Running the APBS calculations, which is the same as ::

          apbs apbs.in

        Parameters
        ----------
        apbs_exe : str, optional
            The path to the APBS program.
            (``/opt/local/bin/apbs``)
        apbs_in: str, optional
            The input file to the APBS program. (``apbs.in``)
        apbs_out : str, optional
            The output file for the APBS calculation. (``apbs.out``)
        '''
        call([apbs_exe, apbs_in])

    def _dx2IP(self, dx):
        g = Grid(dx)
        V = ((np.prod(g.delta*pq.angstrom)) * np.prod(g.grid.shape))
        self.apbs_vol = V
        return np.average(g.grid) * (1/pq.e) * (constants.kB * self.temp) * V

    def read_APBS(self, ligand_RIP_het='ligand_RIP_het.dx',
                  protein_RIP_het='protein_RIP_het.dx',
                  ligand_RIP_hom='ligand_RIP_hom.dx', IP=None):
        ''' Read the result from the APBS calculations.

        Parameters
        ----------
        ligand_RIP_het : str, optional
            The APBS program output (``ligand_RIP_het.dx``).
        protein_RIP_het : str, optional
            The APBS program output (``protein_RIP_het.dx``).
        ligand_RIP_hom : str, optional
            The APBS program output (``ligand_RIP_hom.dx``).
        IP : float, optional
            If system only has ligand, set the integrated potential of protein
            to 0.
        '''
        # Ligand Het
        IL_Bx = self._dx2IP(ligand_RIP_het)
        IL_BQx = (-constants.xi_CB * constants.coulomb_factor / self.water.epsilon_S) * self.apbs_qL * (self.apbs_vol ** (2.0 / 3.0))
        self.IL = IL_Bx - IL_BQx
        # Protein Het
        if IP is not None:
            self.IP = IP
        if self.IP is None:
            IP_Bx = self._dx2IP(protein_RIP_het)
            IP_BQx = (-constants.xi_CB * constants.coulomb_factor / self.water.epsilon_S) * self.apbs_qP * (self.apbs_vol ** (2.0 / 3.0))
            self.IP = IP_Bx - IP_BQx
        else:
            self.IP = pq.Quantity(self.IP, self.IL.units)
        # Ligand Het
        IL_hom_Bx = self._dx2IP(ligand_RIP_hom)
        IL_hom_BQx = (-constants.xi_CB * constants.coulomb_factor / 1) * self.apbs_qL * (self.apbs_vol ** (2.0 / 3.0))
        IL_hom = IL_hom_Bx - IL_hom_BQx
        self.IL_SLV = self.IL - IL_hom


    def compute(self, NS=None):
        ''' Compute the result.

        Parameters
        ----------
        NS : int, optional
            Rest the number of solvent molecules.

        Returns
        -------
        results : float
            The total correction free energy in cal/mol.
        '''
        if NS:
            self.NS = NS
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
        '''Write the decomposed results

        Parameters
        ----------
        outfile : str
            The output file name.'''
        with open(outfile, 'w') as f:
            f.write('\n'.join(self.output))
