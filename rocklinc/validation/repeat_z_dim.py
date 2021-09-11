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
lipid_model.repeats = 7
# for dim in [97, 129, 161, 193, 225, 257]:
for dim in [449, 481, 513, 545, 577]:
    out_ESP, out_charge, out_diel = lipid_model.run_apbs((300, 300, 100),
                                                         (33,  33, dim))
    ESP_list.append(out_ESP)
# np.save('z_dim.npy', np.vstack(ESP_list))
for i, dim in enumerate([449, 481, 513, 545, 577]):
    plt.plot(np.linspace(0,50,dim), ESP_list[i]-min(ESP_list[i]), label=dim)
plt.xlim(20,30)
plt.ylim(0, 0.2)
plt.legend()
plt.show()