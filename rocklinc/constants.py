'''Various constants'''

import quantities as pq

# physical constants
xi_CB = -2.38008 #Coulomb integration constant see eq. (21) in paper and accompanying text
coulomb_factor = 138.93545585 * (1000 * pq.J * pq.nanometer / pq.e ** 2 / pq.mole)  # (4*pi*eps_0)**-1 in (kJ nm)/(e**2 mol)
kB = 0.0083144621 * (1000 * pq.J / pq.mole/ pq.kelvin) # kJ/(mol K)
# xi LS is the cubic lattice-sum (Wigner) integration constant
xi_LS = -2.837297
# epsilon_0 is the permittivity of vacuum.
# 1/(4π*epsilon_0) is 138.93545585 kJ nm e-2 mol-1
epsilon_0 = 5.727657570135483e-04 * pq.mol * pq.e ** 2 / ((1000 * pq.J) * pq.nm)