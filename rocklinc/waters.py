'''Different water models and related constants'''
import quantities as pq

class TIP4P():
    # For tip4p=51: https://pubs.rsc.org/en/content/articlelanding/2011/CP/c1cp22168j#!divAbstract
    epsilon_S = 51
    # For TIP4P =  * hydrogen charge * OH bond length ^2 + MW charge * OD bond length ^2
    gamma_s = 2 * 0.52 * pq.e * (0.09572 * pq.nm) ** 2 + -1.04 * pq.e * (
                0.015 * pq.nm) ** 2

class TIP3P():
    # epsilon_S is the static relative dielectric permittivity of the solvent, for TIP3 is 97.
    epsilon_S = 97
    # For TIP3P, it is 0.00764 e nm2 = 2 * hydrogen charge [0.417*pq.e] * OH bond length (0.09572*pq.nm)**2
    gamma_s = 2 * 0.00764 * pq.e * pq.nm ** 2