import matplotlib.pyplot as plt
import numpy as np
from rocklinc.lipids import POPC
from gridData import Grid
lipid_model = POPC()
out_ESP, out_charge, out_diel = lipid_model.run_apbs((80, 80, 80), (257,  257, 257))

grid = Grid('lipid.dx')
plt.imshow(grid.grid[128,:,:], vmin=-1, vmax=1, cmap='bwr')
plt.show()