import quantities as pq
import numpy as np
from gridData import Grid
from pkg_resources import resource_filename
import shutil
from subprocess import call

class LipdBase():
    # epsilon
    epsilon_headgroup = 80
    epsilon_hydrophobic = 2
    
    # Lipid data
    thickness = np.nan * pq.angstrom
    #headgroup
    magnitude = np.nan * pq.e / pq.angstrom ** 3
    neg_loc = np.nan * pq.angstrom
    neg_width = np.nan * pq.angstrom
    pos_loc = np.nan * pq.angstrom
    pos_width = np.nan * pq.angstrom

    # APBS related
    repeats = 1

    def gauss(self, x, mu, sigma, A):
        return A * np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)

    def _bimodal(self, x, magnitude, neg_loc, neg_width, pos_loc, pos_width):
        return self.gauss(x, neg_loc, neg_width, -1 * magnitude) + \
               self.gauss(x, pos_loc, pos_width, magnitude)


    def _generate_charge(self, z_list):
        charge_list = []
        dim = self.dim[-1]
        thickness = self.thickness
        for z in z_list:
            perodic_z = z % dim
            if perodic_z >= dim / 2:
                real_x = perodic_z - dim / 2
                charge_list.append(self._bimodal(real_x, self.magnitude, self.neg_loc,
                                     self.neg_width, self.pos_loc,
                                     self.pos_width))
            elif perodic_z < dim / 2:
                real_x = dim / 2 - perodic_z
                charge_list.append(self._bimodal(real_x, self.magnitude, self.neg_loc,
                                     self.neg_width, self.pos_loc,
                                     self.pos_width))
            # if perodic_z >= dim / 2:
            #     real_x = dim - perodic_z
            #     charge_list.append(self._bimodal(real_x, self.magnitude, self.neg_loc,
            #                         self.neg_width, self.pos_loc,
            #                         self.pos_width))
            # else:
            #     real_x = perodic_z
            #     charge_list.append(self._bimodal(real_x, self.magnitude, self.neg_loc,
            #                         self.neg_width, self.pos_loc,
            #                         self.pos_width))

        # correct the unit
        charge_list = np.array(charge_list) * charge_list[0].units
        return charge_list

    def _generate_diel(self, z_list):
        diel_list = []
        dim = self.dim[-1]
        thickness = self.thickness
        for z in z_list:
            perodic_z = z % dim
            if perodic_z >= dim / 2 and perodic_z <= dim / 2 + thickness:
                diel_list.append(self.epsilon_hydrophobic)
            elif perodic_z <= dim / 2 and perodic_z >= dim / 2 - thickness:
                diel_list.append(self.epsilon_hydrophobic)
            else:
                diel_list.append(self.epsilon_headgroup)
            # if perodic_z <= thickness:
            #     diel_list.append(self.epsilon_hydrophobic)
            # elif perodic_z >= dim - thickness:
            #     diel_list.append(self.epsilon_hydrophobic)
            # else:
            #     diel_list.append(self.epsilon_headgroup)
        return np.array(diel_list)

    def gen_apbs_files(self, dim, grid):
        self.dim = pq.Quantity(dim, pq.angstrom)
        self.grid = grid
        # The grid will include the first and last point
        self.delta = self.dim / (np.array(self.grid) - 1)

        # generate the one dimensional z axis
        z_list = np.linspace(0 * pq.angstrom, self.dim[-1] * self.repeats, self.grid[-1])

        charge = self._generate_charge(z_list)
        diel = self._generate_diel(z_list)

        # expand to 3D
        diel = np.ones((self.grid[0], self.grid[1], 1)) * diel.reshape(
            (1, 1, self.grid[2]))

        charge = np.ones((self.grid[0], self.grid[1], 1)) * charge.reshape(
            (1, 1, self.grid[2]))
        # Scale unit to e
        charge = charge #* self.delta[0] * self.delta[1] * self.delta[2]

        origin = self.dim / 2 * -1
        origin = origin.rescale(pq.angstrom)
        delta = self.delta.rescale(pq.angstrom)

        # Write diel
        # write x
        origin_x = origin.copy()
        origin_x[0] += self.delta[0] / 2
        g = Grid(diel, origin=origin_x.magnitude, delta=delta.magnitude)
        g.export('dielx.dx', typequote='')

        # write y
        origin_y = origin.copy()
        origin_y[1] += self.delta[1] / 2
        g = Grid(diel, origin=origin_y.magnitude, delta=delta.magnitude)
        g.export('diely.dx', typequote='')

        # write z
        origin_z = origin.copy()
        origin_z[2] += self.delta[2] / 2
        g = Grid(diel, origin=origin_z.magnitude, delta=delta.magnitude)
        g.export('dielz.dx', typequote='')

        # Write charge
        charge = charge.rescale(pq.e/pq.angstrom**3)
        g = Grid(charge.magnitude, origin=origin.magnitude, delta=delta.magnitude)
        g.export('charge.dx', typequote='')

    def gen_apbs_input(self):
        with open(resource_filename(__name__, 'data/lipid.in'), 'r') as f:
            txt = f.read()
        with open('lipid.in', 'w') as f:
            f.write(txt.format(grid_x=self.grid[0],
                               grid_y=self.grid[1],
                               grid_z=self.grid[2],
                               x=self.dim[0].magnitude,
                               y=self.dim[1].magnitude,
                               z=self.dim[2].magnitude,
                               ))
        shutil.copy(resource_filename(__name__, 'data/ref.pqr'), './')

    def run_apbs(self, dim, grid):
        self.gen_apbs_files(dim, grid)
        self.gen_apbs_input()
        call('/opt/local/bin/apbs lipid.in', shell=True)

        charge = Grid('charge_check.dx').grid
        diel = Grid('dielx_check.dx').grid
        lipid = Grid('lipid.dx').grid
        out_ESP = lipid[int(np.rint(self.grid[0]/2)), int(np.rint(self.grid[1]/2)), :]
        out_charge = charge[int(np.rint(self.grid[0]/2)), int(np.rint(self.grid[1]/2)), :]
        out_diel = diel[int(np.rint(self.grid[0]/2)), int(np.rint(self.grid[1]/2)), :]

        return out_ESP, out_charge, out_diel


class POPC(LipdBase):
    # hydrophobic thickness
    thickness = 13.80515 * pq.angstrom
    #headgroup
    magnitude = 6.22602388e-04 * pq.e / pq.angstrom ** 3
    neg_loc = 18.5969598 * pq.angstrom
    neg_width = 1.98288109 * pq.angstrom
    pos_loc = 24.1998191 * pq.angstrom
    pos_width = 1.91347013 * pq.angstrom
