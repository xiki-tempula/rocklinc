from gridData import Grid
import matplotlib.pyplot as plt
import numpy as np
from rocklinc.lipids import POPC

lipid_model = POPC()
out_ESP, out_charge, out_diel = lipid_model.run_apbs((300, 300, 80),
                                                     (257, 257, 257))
plt.plot(np.linspace(0,8,257), out_ESP)
plt.show()
plt.plot(np.linspace(0,8,257), out_charge)
plt.show()
plt.plot(np.linspace(0,8,257), out_diel)
plt.show()