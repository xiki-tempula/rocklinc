from gridData import Grid
import matplotlib.pyplot as plt
import numpy as np
from rocklinc.lipids import POPC
# Obtain the possible dime value
L = 4
# n = c*(2**(L+1)) + 1
res = []
for c in range(20):
    res.append(c * (2 ** (L + 1)) + 1)
# [1, 33, 65, 97, 129, 161, 193, 225, 257, 289, 321, 353, 385, 417, 449, 481, 513, 545, 577, 609]
lipid_model = POPC()
ESP_list = []
for dim in [33, 65, 97, 129, 161, 193, 225, 257]:
    out_ESP, out_charge, out_diel = lipid_model.run_apbs((300, 300, 100),
                                                         (dim,  dim, 257))
    ESP_list.append(out_ESP)
np.save('xy_dim.npy', np.vstack(ESP_list))
for i, dim in enumerate([33, ]):
    plt.plot(np.linspace(0,10,257), ESP_list[i], label=dim)
plt.legend()
plt.show()