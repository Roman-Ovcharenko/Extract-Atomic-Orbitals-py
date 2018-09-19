import numpy as np

from others import get_symb_ll
from others import reconstruct_rad_orbital
from others import integrate

class RedOrbital(object):

    def __init__(self, orbitals, irorb, avg_list, at_symb, bas):
        self.irorb = irorb
        self.avg_list = avg_list
        self.num_deg_orb = len(avg_list)
        self.nn = orbitals[self.avg_list[0]].nn
        self.ll = orbitals[self.avg_list[0]].ll
        self.jj = orbitals[self.avg_list[0]].jj
        self.orb_sym = self._get_orb_symm(at_symb)
        self.energy = self._avg_energy_over_mj(orbitals)
        self.exps = bas.exps[self.ll]
        self.cf = self._avg_cf_rad_mj(orbitals)
        self.Rc = orbitals[self.avg_list[0]].Rc
        self.phi = reconstruct_rad_orbital(self.cf, self.exps, self.Rc, self.ll)
        self._rad_orbitals_printout()
        return None

    def _rad_orbitals_printout(self):
        norm = integrate(self.Rc, self.phi, self.phi)
        if self.irorb == 0:
            print("\n  irorb       orb_sym           norm")
        print(" {:>4d}        {:<14s}   {:6.3f}".format(self.irorb, self.orb_sym, norm)) 
        file_ = self.orb_sym + ".wf"
        with open(file_, "w") as f:
            for ir in range(len(self.Rc)):
                f.write("  {:22.16e}    {:22.16e}\n".format(self.Rc[ir], self.phi[ir]))
        return None

    def _avg_cf_rad_mj(self, orbitals):
        cf_avg = [0.0 for iexp in range(len(self.exps))]
        for iexp in range(len(self.exps)):
            sign = np.sign(orbitals[self.avg_list[0]].cf_rad[iexp])
            for iorb in range(self.num_deg_orb):
                cf_avg[iexp] = (cf_avg[iexp] 
                                + sign * abs(orbitals[self.avg_list[iorb]].cf_rad[iexp]))
            cf_avg[iexp] = cf_avg[iexp] / self.num_deg_orb
        return cf_avg

    def _get_orb_symm(self, at_symb):
        orb_sym = ("{:s}-{:d}{:s}_{:.1f}"
                    .format(at_symb, self.nn, get_symb_ll(self.ll), self.jj))
        return orb_sym

    def _avg_energy_over_mj(self, orbitals):
        energy = 0.0
        for iorb in range(self.num_deg_orb):
            energy = energy + orbitals[self.avg_list[iorb]].energy
        energy = energy / self.num_deg_orb
        return energy



