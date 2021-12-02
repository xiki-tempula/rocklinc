import quantities as pq
import numpy as np
from gridData import Grid
from pkg_resources import resource_filename
import shutil
from subprocess import call
import matplotlib.pyplot as plt

def gauss(x, mu, sigma, A):
    return A * np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)

class LipdBase():
    ''' setup the grids for the APBS calculations.

    Parameters
    ----------
    dim : list
        The unitcell dimensions of the system in Å ``[lx, ly, lz]``.
        The pq.Quantity array could also be used.
    grid: list
        The number of grids in the x, y, z axis ``[257, 257, 257]``.
    '''
    def __init__(self, dim, grid):
        self.dim = pq.Quantity(dim, pq.angstrom)
        self.grid = grid
        # The grid will include the first and last point
        self.delta = self.dim / (np.array(self.grid) - 1)

    @staticmethod
    def gaussian_mixture(x, **kwargs):
        ''' Charge density with respect to the distance to the end of the
        lipid.

        Parameters
        ----------
        x : float or array
            The distance to the end of the lipid tail.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        charge : array
            The one-dimensional array of the charge density.
        '''
        ...

    def generate_charge(self, z_list, param_charge):
        ''' Generate the charge density from the coordinates.

        Parameters
        ----------
        z_list : array
            The one-dimensional coordinate array.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        charge : array
            The one-dimensional array of the charge density.
        '''
        charge_list = []
        dim = self.dim[-1]
        for z in z_list:
            perodic_z = z % dim
            if perodic_z >= dim / 2:
                real_x = perodic_z - dim / 2
                charge_list.append(self.gaussian_mixture(real_x, *param_charge))
            else:
                real_x = dim / 2 - perodic_z
                charge_list.append(self.gaussian_mixture(real_x, *param_charge))

        # correct the unit
        charge_list = np.array(charge_list) * charge_list[0].units
        return charge_list

    @staticmethod
    def diel_model(x, **kwargs):
        ''' The dielectric constant with respect to the distance to the end
        of the lipid.

        Parameters
        ----------
        x : float or array
            The distance to the end of the lipid tail.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        diel : array
            The one-dimensional array of the dielectric constant.
        '''
        ...

    def generate_diel(self, z_list, param_diel):
        ''' Generate the dielectric constant from the coordinates.

        Parameters
        ----------
        z_list : array
            The one-dimensional coordinate array.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        diel : array
            The one-dimensional array of the dielectric constant.
        '''
        diel_list = []
        dim = self.dim[-1]
        for z in z_list:
            perodic_z = z % dim
            if perodic_z >= dim / 2:
                real_x = perodic_z - dim / 2
                diel_list.append(self.diel_model(real_x, *param_diel))
            else:
                real_x = dim / 2 - perodic_z
                diel_list.append(self.diel_model(real_x, *param_diel))
        return np.array(diel_list)

    def write_dx_input(self, param_charge, param_diel):
        ''' Generate the dielectric constant and charge grid for APBS
        calculations.

        ``generate_charge`` and ``generate_diel`` are used to generate the
        profile of the charge and dielectric constant on the one-dimensional
        z-axis which is then broadcast to cover the whole 3D space. The
        generated files are named as ``dielx.dx``, ``diely.dx``,
        ``dielz.dx`` and ``charge.dx``.

        Parameters
        ----------
        **kwargs : dict
            ``generate_charge`` uses the arguments from kwargs['charge'] and
            ``generate_diel`` uses the arguments from kwargs['diel'].
        '''
        # generate the one dimensional z axis
        z_list = np.linspace(0 * pq.angstrom, self.dim[-1], self.grid[-1])

        charge = self.generate_charge(z_list, param_charge)
        diel = self.generate_diel(z_list, param_diel)

        # expand to 3D
        diel = np.ones((self.grid[0], self.grid[1], 1)) * diel.reshape(
            (1, 1, self.grid[2]))

        charge = np.ones((self.grid[0], self.grid[1], 1)) * charge.reshape(
            (1, 1, self.grid[2]))

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

    def write_apbs_input(self):
        ''' Generate the APBS input file ``lipid.in`` and the reference qpr
        file ``ref.pqr``.

        The APBS calculation will use the ``dielx.dx``, ``diely.dx`` and
        ``dielz.dx`` for dielectric constant and ``charge.dx`` for the
        charge density. The ``ref.pqr`` is a place holder and is not used in
        the calculation.
        '''
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

    def run_APBS(self, apbs_exe='/opt/local/bin/apbs', apbs_in='lipid.in'):
        ''' Running the APBS calculations, which is the same as ::

          apbs lipid.in

        Parameters
        ----------
        apbs_exe : str, optional
            The path to the APBS program.
            (``/opt/local/bin/apbs``)
        apbs_in: str, optional
            The input file to the APBS program. (``lipid.in``)
        '''
        call([apbs_exe, apbs_in])

    def read_APBS(self):
        ''' Read the APBS output.

        Returns
        -------
        charge : array
            The three-dimensional array of the charge density.
        diel : array
            The three-dimensional array of the dielectric constant in the x
            direction.
        lipid : array
            The three-dimensional array of the ESP.
        '''
        charge = Grid('charge_check.dx').grid
        diel = Grid('dielx_check.dx').grid
        lipid = Grid('lipid.dx').grid
        return charge, diel, lipid

    def read_central_ESP(self):
        ''' Read the one-dimensional ESP along the Z-axis at the center of
        the x-y plane.

        Returns
        -------
        out_ESP : array
            The one-dimensional ESP.
        '''
        lipid = Grid('lipid.dx').grid
        out_ESP = lipid[int(np.rint(self.grid[0] / 2)),
                  int(np.rint(self.grid[1] / 2)), :]
        return out_ESP

class StepLipid(LipdBase):
    @staticmethod
    def gaussian_mixture(x, magnitude, neg_loc, neg_width, pos_loc, pos_width):
        ''' Charge density with respect to the distance to the end of the
        lipid.

        This functional form is the sum of two Gaussians with opposite
        magnitude.

        Parameters
        ----------
        x : float or array
            The distance to the end of the lipid tail.
        magnitude : float
            The absolute magnitude of both Gaussians (e/Å^3).
        neg_loc : float
            The position of the Gaussians with negative magnitude (Å).
        neg_width : float
            The spread of the Gaussians with negative magnitude (Å).
        pos_loc : float
            The position of the Gaussians with positive magnitude (Å).
        pos_width : float
            The spread of the Gaussians with positive magnitude (Å).

        Returns
        -------
        charge : array
            The one-dimensional array of the charge density.
        '''
        return gauss(x, pq.Quantity(neg_loc, pq.angstrom),
                     pq.Quantity(neg_width, pq.angstrom),
                     -1 * pq.Quantity(magnitude, pq.e / pq.angstrom ** 3)) + \
               gauss(x, pq.Quantity(pos_loc, pq.angstrom),
                     pq.Quantity(pos_width, pq.angstrom),
                     pq.Quantity(magnitude, pq.e / pq.angstrom ** 3))

    @staticmethod
    def diel_model(x, thickness, diel_low, diel_high):
        ''' The dielectric constant with respect to the distance to the end
        of the lipid.

        Parameters
        ----------
        x : float or array
            The distance to the end of the lipid tail.
        thickness : float
            The thickness of the lipid acryl chain, which is the distance from
            the end of the acryl chain and the begining of the acryl chain.
        diel_low : float
            The dielectric constant of the region of the lipid acryl chain.
        diel_high : float
            The dielectric constant of the region of the lipid head group
            and solvent.

        Returns
        -------
        diel : array
            The one-dimensional array of the dielectric constant.
        '''
        if x < thickness:
            return diel_low
        else:
            return diel_high

class CurveLipid(LipdBase):
    @staticmethod
    def gaussian_mixture(x, cho_loc, cho_std, cho_mag,
                         po4_loc, po4_std, po4_mag,
                         GL_loc, GL_std, GL_mag,
                         neg_loc, neg_std, neg_mag,
                         pos_std, pos_mag):
        ''' Charge density with respect to the distance to the end of the
        lipid.

        This functional form is the sum of two Gaussians with opposite
        magnitude.

        Parameters
        ----------
        x : float or array
            The distance to the end of the lipid tail.
        cho_loc : float
            The position of the Gaussians of the choline group (Å).
        cho_std : float
            The spread of the Gaussians of the choline group (Å).
        cho_mag : float
            The absolute magnitude of choline group (e/Å^3).
        po4_loc : float
            The position of the Gaussians of the phosphate group (Å).
        po4_std : float
            The spread of the Gaussians of the phosphate group (Å).
        po4_mag : float
            The absolute magnitude of phosphate group (e/Å^3).
        GL_loc : float
            The position of the Gaussians of the phosphate group (Å).
        GL_std : float
            The spread of the Gaussians of the phosphate group (Å).
        GL_mag : float
            The absolute magnitude of phosphate group (e/Å^3).
        neg_loc : float
            The position of the Gaussians of the negative dipole of acryl chain (Å).
        neg_std : float
            The spread of the Gaussians of negative dipole of acryl chain (Å).
        neg_mag : float
            The absolute magnitude of negative dipole of acryl chain (e/Å^3).
        pos_std : float
            The spread of the Gaussians with positive dipole of acryl chain (Å).
        pos_mag : float
            The absolute magnitude of positive dipole of acryl chain (e/Å^3).

        Returns
        -------
        charge : array
            The one-dimensional array of the charge density.
        '''
        return gauss(x, pq.Quantity(cho_loc, pq.angstrom),
                     pq.Quantity(cho_std, pq.angstrom),
                     pq.Quantity(cho_mag, pq.e / pq.angstrom ** 3)) + \
               gauss(x, pq.Quantity(po4_loc, pq.angstrom),
                     pq.Quantity(po4_std, pq.angstrom),
                     pq.Quantity(po4_mag, pq.e / pq.angstrom ** 3) * -1) + \
               gauss(x, pq.Quantity(GL_loc, pq.angstrom),
                     pq.Quantity(GL_std, pq.angstrom),
                     pq.Quantity(GL_mag, pq.e / pq.angstrom ** 3)) + \
               gauss(x, pq.Quantity(neg_loc, pq.angstrom),
                     pq.Quantity(neg_std, pq.angstrom),
                     pq.Quantity(neg_mag, pq.e / pq.angstrom ** 3) * -1) + \
               gauss(x, pq.Quantity(0, pq.angstrom),
                     pq.Quantity(pos_std, pq.angstrom),
                     pq.Quantity(pos_mag, pq.e / pq.angstrom ** 3))

    @staticmethod
    def diel_model(x, k1, b1, k2, b2):
        ''' The dielectric constant with respect to the distance to the end
        of the lipid.

        Parameters
        ----------
        x : float or array
            The distance to the end of the lipid tail.
        k1 : float
            The slope of the switching from the low to the high dielectric
            constant region. (1/Å)
        b1 : float
            The length of the acryl chain. (Å)
        k2 : float
            The dielectric constant of the solvent minus the dielectric
            constant of the acryl chain of the lipid.
        b2 : float
            The dielectric constant of the acryl chain of the lipid.

        Returns
        -------
        diel : array
            The one-dimensional array of the dielectric constant.
        '''
        return k2 / (1 + np.exp(-(pq.Quantity(k1, 1/pq.angstrom)*(x-pq.Quantity(b1, pq.angstrom))))) + b2

def plot_dx(x, charge, diel, esp):
    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()

    # Offset the right spine of par2.  The ticks and label have already been
    # placed on the right by twinx above.
    par2.spines["right"].set_position(("axes", 1.2))
    # Having been created by twinx, par2 has its frame off, so the line of its
    # detached spine is invisible.  First, activate the frame but make the patch
    # and spines invisible.
    par2.set_frame_on(True)
    par2.patch.set_visible(False)
    for sp in par2.spines.values():
        sp.set_visible(False)
    # Second, show the right spine.
    par2.spines["right"].set_visible(True)

    shape = charge.shape
    charge = charge[int(np.rint(shape[0] / 2)), int(np.rint(shape[0] / 2)), :]
    diel = diel[int(np.rint(shape[0] / 2)), int(np.rint(shape[0] / 2)), :]
    esp = esp[int(np.rint(shape[0] / 2)), int(np.rint(shape[0] / 2)), :]

    p1, = host.plot(x, esp, "b-", label="ESP")
    p2, = par1.plot(x, charge*1000, "r-", label="Charge")
    p3, = par2.plot(x, diel, "g-", label="Diel")

    host.set_xlabel("Z axis (A)")
    host.set_ylabel("ESP (kT/e)")
    par1.set_ylabel("Charge (me/Å^3)")
    par2.set_ylabel("Diel")

    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    par2.yaxis.label.set_color(p3.get_color())

    tkw = dict(size=4, width=1.5)
    host.tick_params(axis='y', colors=p1.get_color(), **tkw)
    par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    host.tick_params(axis='x', **tkw)

    lines = [p1, p2, p3]

    host.legend(lines, [l.get_label() for l in lines])
    return fig