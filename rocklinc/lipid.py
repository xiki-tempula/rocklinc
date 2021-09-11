import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from gridData import Grid
from subprocess import call

diel_membrane = 2
diel_sol = 80
class NewLipid():
    def __init__(self, universe, lipid_selection='resname POPC'):
        self.universe = universe
        self.lipid = universe.select_atoms(lipid_selection)



    def _bimodal(self, x, A, mu1, sigma1, mu2, sigma2):
        def gauss(x, mu, sigma, A):
            return A * np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)

        # Set the amplitude of the two Gaussians to be the same
        # Such that it will have a neutral charge
        return gauss(x, mu1, sigma1, -1 * A) + gauss(x, mu2, sigma2, A)

    def _fit_bimodal_gauss(self, x, y, guess):
        params, cov = curve_fit(self._bimodal, x, y, guess)
        return params

    def cal_hydrophobic_thickness(self, name='name C32 C22'):
        thickness_list = []
        for ts in self.universe.trajectory:
            z = self.lipid.select_atoms(name).positions[:, -1]
            middle = np.mean(z)
            thickness = np.mean(np.abs(z-middle))
            thickness_list.append(thickness)
        self.thickness = thickness_list

    def generate_model(self, initial_delta=0.4, guess=None):
        z_list = []
        box = []
        # Center the center of the hydrophobic part to z=0
        for ts in self.universe.trajectory:
            self.lipid.translate([0, 0, -self.lipid.center_of_mass()[-1]])
            z = np.abs(self.lipid.positions[:, -1])
            z_list.append(z)
            box.append(ts.dimensions[:2])

        # derive the xy density
        xy_plane = np.mean(box, axis=0)
        vol = xy_plane[0] * xy_plane[1]

        # Initial fitting
        # set up the inital grid
        max_length = max(np.hstack(z_list))
        histogram = [[] for i in range(int(max_length / initial_delta) + 1)]

        charges = self.lipid.charges
        for frame in z_list:
            for index in range(len(charges)):
                histogram[int(frame[index] / initial_delta)].append(charges[index])

        num_frame = len(self.universe.trajectory)

        sum_charge = np.array([np.sum(line) for line in histogram])
        # Pre frame, pre A^2 in the xy-plane, pre A in the z-axis, per leaflet
        unit_charge = sum_charge / num_frame / vol / initial_delta / 2

        x = np.arange(len(unit_charge)) * initial_delta
        y = unit_charge

        if guess is None:
            guess = [0.0006, 18.8, 2.0, 24.3, 2.0]

        initial_fit = self._fit_bimodal_gauss(x, y, guess)
        _, _, std1, _, std2 = initial_fit
        # Set the optimal delta to 1/5 th of the spread
        delta = min(std1 / 5, std2 / 5)

        # second fit
        histogram = [[] for i in range(int(max_length / delta) + 1)]
        for frame in z_list:
            for index in range(len(charges)):
                histogram[int(frame[index] / delta)].append(charges[index])

        sum_charge = np.array([np.sum(line) for line in histogram])

        unit_charge = sum_charge / num_frame / vol / delta / 2

        x = np.arange(len(unit_charge)) * delta
        y = unit_charge

        self.params = self._fit_bimodal_gauss(x, y, initial_fit)
        self.point_charge_range = x
        self.point_charge_hist = y
        return self.params

    def plot_model(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(self.point_charge_range, self.point_charge_hist)
        x = np.linspace(min(self.point_charge_range),
                        max(self.point_charge_range),
                        100)
        ax.plot(x, self._bimodal(x,*self.params),color='red',lw=3,label='model')
        ax.set_xlabel('Distance (Å)')
        ax.set_ylabel('Charge density (e/Å^3)')
        ax.set_title(
            'Amplitude: {:.4f}; peak 1 {:.1f}±{:.1f} peak 2 {:.1f}±{:.1f}'.format(
                *self.params
            ))

        try:
            ax.fill_betweenx([-self.params[0], self.params[0]],
                             np.mean(self.thickness) + np.std(self.thickness),
                             np.mean(self.thickness) - np.std(self.thickness),
                             color='#D2B9D3', zorder=-1)
        except:
            pass
        return ax

    def generate_APBS_input(self, dime=None, grid=None, thickness=None):
        if dime is None:
            dime = [100, ] * 3

        if grid is None:
            grid = (257, 257, 257)

        delta = np.array(dime) / (np.array(grid)-1)
        if thickness is None:
            thickness = self.thickness

        diel, charge = self._APBS_mem_water_mem(dime, grid, thickness)

        # write the dx file
        limit = np.array(dime) / -2

        # write x
        origin_x = limit.copy()
        origin_x[0] += delta[0] / 2
        g = Grid(diel, origin=origin_x, delta=delta)
        g.export('dielx.dx', typequote='')

        # write y
        origin_y = limit.copy()
        origin_y[1] += delta[1] / 2
        g = Grid(diel, origin=origin_y, delta=delta)
        g.export('diely.dx', typequote='')

        # write z
        origin_z = limit.copy()
        origin_z[2] += delta[2] / 2
        g = Grid(diel, origin=origin_z, delta=delta)
        g.export('dielz.dx', typequote='')

        # write charge
        charge = np.ones((grid[0], grid[1], 1)) * charge.reshape((1,1,len(charge)))
        g = Grid(charge, origin=limit, delta=delta)
        g.export('charge.dx', typequote='')

    def _APBS_water_mem_water(self, dime, grid, thickness):
        delta = np.array(dime) / (np.array(grid) - 1)
        # the lipid is perpendicular to the xy plane.
        num_thickness = np.mean(thickness) * 2 / delta[-1]
        start, end = grid[-1] // 2 - num_thickness // 2, grid[
            -1] // 2 + num_thickness // 2

        diel = np.ones(grid) * diel_sol
        diel[:, :, int(start): int(end)] = diel_membrane

        # make x axis
        x = np.linspace(0, dime[-1], grid[-1]) - dime[-1]/2
        A, mu1, sigma1, mu2, sigma2 = self.params
        # unit is e/A^3
        charge = self._bimodal(x, A, mu1, sigma1, mu2, sigma2) + \
                 self._bimodal(x, A, -1*mu1, sigma1, -1*mu2, sigma2)

        return diel, charge

    def _APBS_mem_water_mem(self, dime, grid, thickness):
        delta = np.array(dime) / (np.array(grid) - 1)
        # the lipid is perpendicular to the xy plane.
        num_thickness = int(np.mean(thickness) / delta[-1])

        diel = np.ones(grid) * 80
        diel[:, :, :num_thickness] = 2
        diel[:, :, -num_thickness:] = 2

        # make x axis
        x = np.linspace(0, dime[-1], grid[-1])
        A, mu1, sigma1, mu2, sigma2 = self.params
        # unit is e/A^3
        forward = self._bimodal(x, A, mu1, sigma1, mu2, sigma2)
        charge = forward + np.flip(forward)
        return diel, charge

    def benchmark_xy(self):
        total_list = []
        for xy in range(100, 300, 20):
            self.generate_APBS_input(dime=(xy, xy, 100))
            with open('apbs_templete.in', 'r') as r:
                with open('apbs.in', 'w') as f:
                    f.write(r.read().format(x=xy, y=xy, z=100))
            call('/opt/local/bin/apbs apbs.in', shell=True)
            g = Grid('lipid.dx')
            total_list.append(g.grid[129,129,:])
        return total_list

    def benchmark_z(self):
        total_list = []
        for z in range(100, 300, 20):
            self.generate_APBS_input(dime=(300, 300, z))
            with open('apbs_templete.in', 'r') as r:
                with open('apbs.in', 'w') as f:
                    f.write(r.read().format(x=300, y=300, z=z))
            call('/opt/local/bin/apbs apbs.in', shell=True)
            g = Grid('lipid.dx')
            total_list.append(g.grid[129,129,:])
        return total_list

    def benchmark_sandwich_z(self):
        esp_list = []
        diel_list = []
        charge_list = []
        for z in range(100, 300, 20):
            self.generate_APBS_input(dime=(300, 300, z))
            with open('apbs_templete.in', 'r') as r:
                with open('apbs.in', 'w') as f:
                    f.write(r.read().format(x=300, y=300, z=z))
            call('/opt/local/bin/apbs apbs.in', shell=True)
            g = Grid('lipid.dx')
            esp_list.append(g.grid[129,129,:])
            g = Grid('dielx_check.dx')
            diel_list.append(g.grid[129,129,:])
            g = Grid('charge_check.dx')
            charge_list.append(g.grid[129,129,:])
        return esp_list, diel_list, charge_list

    def benchmark_layer(self):
        # separate = 10 nm
        dime = np.array((300, 300, 100))
        grid = np.array((257, 257, 257))
        thickness = np.mean(self.thickness)

        def diel_value(x):
            perodic_x = x % dime[-1]
            if perodic_x >= dime[-1] / 2 and perodic_x <= dime[-1] / 2 + thickness:
                return diel_membrane
            elif perodic_x <= dime[-1] / 2 and perodic_x >= dime[-1] / 2 - thickness:
                return diel_membrane
            else:
                return diel_sol

        def charge_value(x):
            perodic_x = x % dime[-1]
            if perodic_x >= dime[-1] / 2:
                real_x = perodic_x - dime[-1] / 2
                return self._bimodal(real_x, *self.params)
            elif perodic_x < dime[-1] / 2:
                real_x = dime[-1] / 2 - perodic_x
                return self._bimodal(real_x, *self.params)

        esp_list = []
        diel_list = []
        charge_list = []

        for r in range(1,10):
            x_list = np.linspace(0, dime[-1]*r, grid[-1])
            diel = [diel_value(x) for x in x_list]
            charge = [charge_value(x) for x in x_list]

            current_dim = np.array([dime[0], dime[1], dime[2]*r])
            delta = current_dim/(grid-1)

            fig, ax = plt.subplots(nrows=2)
            ax[0].plot(diel)
            ax[1].plot(charge)
            plt.show()

            diel = np.ones((grid[0], grid[1], 1)) * np.array(diel).reshape(
                (1, 1, grid[2]))
            charge = np.ones((grid[0], grid[1], 1)) * np.array(charge).reshape(
                (1, 1, grid[2]))

            limit = current_dim/2*-1
            # write x
            origin_x = limit.copy()
            origin_x[0] += delta[0] / 2
            g = Grid(diel, origin=origin_x, delta=delta)
            g.export('dielx.dx', typequote='')

            # write y
            origin_y = limit.copy()
            origin_y[1] += delta[1] / 2
            g = Grid(diel, origin=origin_y, delta=delta)
            g.export('diely.dx', typequote='')

            # write z
            origin_z = limit.copy()
            origin_z[2] += delta[2] / 2
            g = Grid(diel, origin=origin_z, delta=delta)
            g.export('dielz.dx', typequote='')

            # write charge
            g = Grid(charge, origin=limit, delta=delta)
            g.export('charge.dx', typequote='')

            with open('apbs_templete.in', 'r') as r:
                with open('apbs.in', 'w') as f:
                    f.write(r.read().format(x=current_dim[0],
                                            y=current_dim[1],
                                            z=current_dim[2]))
            call('/opt/local/bin/apbs apbs.in', shell=True)
            g = Grid('lipid.dx')
            esp_list.append(g.grid[129,129,:])
            g = Grid('dielx_check.dx')
            diel_list.append(g.grid[129,129,:])
            g = Grid('charge_check.dx')
            charge_list.append(g.grid[129,129,:])

        return esp_list, diel_list, charge_list


    def check_water_thickness(self):
        grid = np.array((257, 257, 257))
        thickness = np.mean(self.thickness)
        def diel_value(x):
            perodic_x = x % dime[-1]
            if perodic_x >= dime[-1] / 2 and perodic_x <= dime[-1] / 2 + thickness:
                return diel_membrane
            elif perodic_x <= dime[-1] / 2 and perodic_x >= dime[-1] / 2 - thickness:
                return diel_membrane
            else:
                return diel_sol

        def charge_value(x):
            perodic_x = x % dime[-1]
            if perodic_x >= dime[-1] / 2:
                real_x = perodic_x - dime[-1] / 2
                return self._bimodal(real_x, *self.params)
            elif perodic_x < dime[-1] / 2:
                real_x = dime[-1] / 2 - perodic_x
                return self._bimodal(real_x, *self.params)

        charge_list = []

        for water in [80, 85, 90, 95, 100, 200, 300]:
            dime = np.array((300, 300, water))
            r=3
            x_list = np.linspace(0, dime[-1]*r, grid[-1])
            diel = [diel_value(x) for x in x_list]
            charge = [charge_value(x) for x in x_list]

            current_dim = np.array([dime[0], dime[1], dime[2]*r])
            delta = current_dim/(grid-1)

            fig, ax = plt.subplots(nrows=2)
            ax[0].plot(diel)
            ax[1].plot(charge)
            plt.show()

            diel = np.ones((grid[0], grid[1], 1)) * np.array(diel).reshape(
                (1, 1, grid[2]))
            charge = np.ones((grid[0], grid[1], 1)) * np.array(charge).reshape(
                (1, 1, grid[2]))

            limit = current_dim/2*-1
            # write x
            origin_x = limit.copy()
            origin_x[0] += delta[0] / 2
            g = Grid(diel, origin=origin_x, delta=delta)
            g.export('dielx.dx', typequote='')

            # write y
            origin_y = limit.copy()
            origin_y[1] += delta[1] / 2
            g = Grid(diel, origin=origin_y, delta=delta)
            g.export('diely.dx', typequote='')

            # write z
            origin_z = limit.copy()
            origin_z[2] += delta[2] / 2
            g = Grid(diel, origin=origin_z, delta=delta)
            g.export('dielz.dx', typequote='')

            # write charge
            g = Grid(charge, origin=limit, delta=delta)
            g.export('charge.dx', typequote='')

            with open('apbs_templete.in', 'r') as r:
                with open('apbs.in', 'w') as f:
                    f.write(r.read().format(x=current_dim[0],
                                            y=current_dim[1],
                                            z=current_dim[2]))
            call('/opt/local/bin/apbs apbs.in', shell=True)
            g = Grid('lipid.dx')
            charge_list.append(g.grid[129,129,:])
        return charge_list
