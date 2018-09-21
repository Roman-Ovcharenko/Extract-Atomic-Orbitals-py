from math import factorial as fact
from math import sqrt
from math import log
from math import exp
from math import pi
from numpy import sign

from others import binom
from others import clebsch
from others import integrate
from others import get_symb_ll
from others import reconstruct_rad_orbital
from others import dfactorial as dfact

##################################################################################################
# The Orbital object represents an atomic orbital
#
# self.file_prop - the file containing coefficients of atomic orbitals
# self.energy - energy of orbital
# self.parity - parity of orbital
# self.iorb - serial number of the given orbital in atom
# self.nn, self.ll, self.jj, self.mj - quantum numbers 
# self.cf_La, self.cf_Lb - orbital coefficients for spin up = a and spin down = b components
# self.orb_sym - symbolic representation of the given orbital
# self.exps - basis set exponents for a given orbital symmetry of the atomic orbital
# self.cf_rad - radial atomic orbital coefficients averaged over spin-up and spin-down 
# configurations
# self.Rc - radial grid for representation of the radial atomic orbital in direct space
# self.phi - radial atomic orbitals represented on the radial grid
# self.norm_avg - normalization of the atomic orbitals averaged over spin-up and spin-down
# configurations
##################################################################################################
class Orbital(object):

    def __init__(self, file_prop, energy, parity, jj, mj, ll, bas, iorb, num_orb_ger, at_symb,
                cf_La, cf_Lb, nn, wrong_norm_orbitals):
        self.file_prop = file_prop
        self.energy = energy
        self.parity = parity
        self.iorb = iorb
        self.mj = mj
        self.jj = jj
        self.cf_La = cf_La
        self.cf_Lb = cf_Lb
        self.ll  = ll
        self.nn = nn
        self.orb_sym = self._get_orb_symm(at_symb)
        self.exps = bas.exps[self.ll]
        cf_rad_up = self._calc_cf_rad(bas, 0.5)
        cf_rad_down = self._calc_cf_rad(bas, -0.5)
        self.cf_rad = self._avg_cf_rad(cf_rad_up, cf_rad_down, bas)
        self.Rc = self._get_radial_grid()
        phi_up = reconstruct_rad_orbital(cf_rad_up, self.exps, self.Rc, self.ll)
        phi_down = reconstruct_rad_orbital(cf_rad_down, self.exps, self.Rc, self.ll)
        self.phi = reconstruct_rad_orbital(self.cf_rad, self.exps, self.Rc, self.ll)
        self.norm_avg = self._rad_orbitals_printout(phi_up, phi_down, self.phi, wrong_norm_orbitals)
        return None

##################################################################################################
# Generate a radial grid for a representation of an atomic orbital in the direct space
##################################################################################################
    def _get_radial_grid(self):
        Rmax = float(20) # AA
        R0 = float(0.0000264588624)  # AA
        NrC = int(15000) 
        dt = log( 1.0 + Rmax/R0 ) / ( NrC - 1.0 )
        Rc = []
        for ir in range(NrC):
            t = ir * dt
            Rc.append(R0 * ( exp(t) - 1.0 ))
        return Rc

##################################################################################################
# Calculation of the coefficient use in the representation of the Cartesian basis function 
# as a spherical harmonic
##################################################################################################
    def _coef(self, lc, lx, ly, lz, gamma):
        if lx + ly + lz != lc:
            raise Exception("lc != lx + ly + lz")
        mult = 2**(lc+9.0/4.0) * gamma**((2*lc+3.0)/4.0) / pi**(1.0/4.0)
        numer = fact(lc+1) * fact(2*lx) * fact(2*ly) * fact(2*lz)
        denum = fact(lx) * fact(ly) * fact(lz) * fact(2*lc+2) 
        return mult * sqrt( numer / denum )

##################################################################################################
# Calculation of the coefficient between Cartesian and spherical Gaussian basis sets 
##################################################################################################
    def _sphr_to_dec_cf(self, lc, ml, lx, ly, lz):
        if lx + ly + lz != lc:
            raise Exception("lc != lx + ly + lz")
        res = complex(0.0, 0.0)
        for lxx in range(lc+1):
            for lyy in range(lc+1):
                for lzz in range(lc+1):
                    if lxx + lyy + lzz == lc:
                        res = res + (self._overlap(lx, ly, lz, lxx, lyy, lzz) 
                                     * self._dec_to_sphr_cf(lc, ml, lxx, lyy, lzz).conjugate())
        return res

##################################################################################################
# Calculation of the overlap matrix
##################################################################################################
    def _overlap(self, lx, ly, lz, lxx, lyy, lzz):
        rlxxx = (lx + lxx) / 2.0
        rlyyy = (ly + lyy) / 2.0
        rlzzz = (lz + lzz) / 2.0
        if rlxxx.is_integer() and rlyyy.is_integer() and rlzzz.is_integer():
            lxxx = int(rlxxx)
            lyyy = int(rlyyy)
            lzzz = int(rlzzz)
        else:
            return 0.0
        numer1 = fact(2*lxxx) * fact(2*lyyy) * fact(2*lzzz)
        denum1 = fact(lxxx) * fact(lyyy) * fact(lzzz)
        numer2 = fact(lx) * fact(ly) * fact(lz) * fact(lxx) * fact(lyy) * fact(lzz) 
        denum2 = fact(2*lx) * fact(2*ly) * fact(2*lz) * fact(2*lxx) * fact(2*lyy) * fact(2*lzz) 
        return numer1 * sqrt(numer2 / denum2) / denum1

##################################################################################################
# Calculation of the coefficient between Cartesian and spherical Gaussian basis sets 
##################################################################################################
    def _dec_to_sphr_cf(self, lc, ml, lx, ly, lz):
        if lx + ly + lz != lc:
            raise Exception("lc != lx + ly + lz")
        mm = abs(ml)
        if ml >= 0:
            sgn_ml = +1
        else:
            sgn_ml = -1
        rj = (lx + ly - mm) / 2.0
        if rj.is_integer():
            j = int(rj)
        else:
            return complex(0.0, 0.0)
        numer = fact(2*lx) * fact(2*ly) * fact(2*lz) * fact(lc) * fact(lc - mm) 
        denum = fact(2*lc) * fact(lx) * fact(ly) * fact(lz) * fact(lc + mm)
        mult1 = sqrt(numer / denum) / (2.0**lc * fact(lc))
        if lc == mm:
            up = 0
        else:
            up = (lc - mm) // 2
        sum_out = complex(0.0, 0.0)
        for ii in range(up+1):
            mult2 = ( binom(lc, ii) * binom(ii, j) * (-1.0)**ii 
                     * fact(2*lc - 2*ii) / fact(lc - mm - 2*ii) ) 
            sum_in = complex(0.0, 0.0)
            for kk in range(j+1):
                sum_in = (sum_in + binom(j, kk) * binom(mm, lx - 2*kk) 
                          * (-1)**(sgn_ml*(mm - lx + 2*kk)/2.0))
            sum_out = sum_out + mult2 * sum_in
        res = mult1 * sum_out
        return res

##################################################################################################
# Construct a symbolic representation for the given atomic orbital
##################################################################################################
    def _get_orb_symm(self, at_symb):
        orb_sym = ("{:s}-{:d}{:s}_{:.1f}_{:+.1f}"
                    .format(at_symb, self.nn, get_symb_ll(self.ll), self.jj, self.mj))
        return orb_sym

##################################################################################################
# Reduce the series of 3D-Cartesian basis set representing atomic orbital into the series of 
# the 1D-Gaussian basis set representing the radial part of atomic orbital 
##################################################################################################
    def _calc_cf_rad(self, bas, beta):
        _eps = 1.0e-10
        exps = bas.exps[self.ll]
        ml = int(self.mj-beta)
        if abs(ml) > self.ll:
            return [ 0.0 for jexp in range(len(exps))]

        if abs(beta - 0.5) < _eps:
            cf = self.cf_La
        elif abs(beta + 0.5) < _eps:
            cf = self.cf_Lb
        else:
            raise Exception("Spin is larger than 0.5")

        sum_ = [complex(0.0, 0.0) for jexp in range(len(exps))]
        for ibas in range(bas.size):
            if bas.ltot[ibas] != self.ll:
                continue
            iexp = exps.index(bas.exp[ibas])
            sum_[iexp] = (sum_[iexp] + cf[ibas] * self._coef(self.ll, bas.lx[ibas], bas.ly[ibas], 
                                                       bas.lz[ibas], bas.exp[ibas])
                          * self._sphr_to_dec_cf(self.ll, ml, bas.lx[ibas], 
                                                 bas.ly[ibas], bas.lz[ibas])) 

        res = [ 0.0 for jexp in range(len(exps))]
        for jexp in range(len(exps)):
            if self.parity == "gerade":
                res[jexp] = sum_[jexp].real 
            else:
                res[jexp] = sum_[jexp].imag 

        return res

##################################################################################################
# Average radial parts for the {n, l, j} orbitals over mj quantum numbers 
##################################################################################################
    def _avg_cf_rad(self, cf_rad_up, cf_rad_down, bas):
        _eps = 1.0e-10
        exps = bas.exps[self.ll]
        beta = 0.5
        cl_up = clebsch(self.jj, self.mj, self.ll, int(self.mj-beta), abs(beta), beta)
        beta = -0.5
        cl_down = clebsch(self.jj, self.mj, self.ll, int(self.mj-beta), abs(beta), beta)
        if (abs(cl_down) < _eps) and (abs(cl_up - 1.0) < _eps):
            cf_rad = cf_rad_up
        elif (abs(cl_down - 1.0) < _eps) and (abs(cl_up) < _eps):
            cf_rad = cf_rad_down
        elif (abs(cl_down) < _eps) and (abs(cl_up) < _eps):
            raise Exception("Wrong quantum numbers:\njj = {}\nmj = {}\nll = {}"
                           .format(self.jj, self.mj, self.ll))
        else:
            cf_rad = []
            for iexp in range(len(exps)):
                cf_rad.append(sign(cf_rad_up[iexp]) * 0.5 * (abs(cf_rad_up[iexp] / cl_up) 
                                                             + abs(cf_rad_down[iexp] / cl_down)))
        return cf_rad

##################################################################################################
# Print the detailed information about radial parts of atomic orbitals to the screen and 
# the radial parts of atomic orbitals themselvev to the corresponding files
##################################################################################################
    def _rad_orbitals_printout(self, phi_up, phi_down, phi_avg, wrong_norm_orbitals):
        norm_up = integrate(self.Rc, phi_up, phi_up)
        norm_down = integrate(self.Rc, phi_down, phi_down)
        norm_avg = integrate(self.Rc, phi_avg, phi_avg)

        beta = 0.5
        cl_up = clebsch(self.jj, self.mj, self.ll, int(self.mj-beta), abs(beta), beta)
        beta = -0.5
        cl_down = clebsch(self.jj, self.mj, self.ll, int(self.mj-beta), abs(beta), beta)

        if ((abs(norm_up - cl_up**2) > 0.1*cl_up**2+0.01) 
            or (abs(norm_down - cl_down**2) > 0.1*cl_down**2+0.01)):
            wrong_norm_orbitals.append(self.orb_sym)

        if self.iorb == 0:
            print("  iorb       orb_sym       norm-up   clebsch-up^2   norm-down"
                  "    clebsch-down^2    norm-avg")
        print(" {:>4d}     {:<14s}   {:6.3f}      {:6.3f}        {:6.3f}        {:6.3f}"
              "           {:6.3f}"
              .format(self.iorb, self.orb_sym, norm_up, cl_up**2, norm_down, cl_down**2,
              norm_avg)) 

        file_ = self.orb_sym + ".wf"

        with open(file_+"-avg", "w") as f_avg, open(file_+"-down", "w") as f_down, open(
            file_+"-up", "w") as f_up:
            for ir in range(len(self.Rc)):
                f_down.write("  {:22.16e}    {:22.16e}\n".format(self.Rc[ir], phi_down[ir]))
                f_up.write("  {:22.16e}    {:22.16e}\n".format(self.Rc[ir], phi_up[ir]))
                f_avg.write("  {:22.16e}    {:22.16e}\n".format(self.Rc[ir], phi_avg[ir]))
        return norm_avg

##################################################################################################
# Calculation of the overlap matrix of two Cartesian Gaussian basis set functions 
##################################################################################################
    def _overlap_dec_basis(self, lxx, lyy, lzz, mll, gamma2, lx, ly, lz, ml, gamma):
        lsum = lx + ly + lz
        lsum2 = lxx + lyy + lzz
        if (lsum != lsum2) or (mll != ml):
            return 0.0
        rad_int = (dfact(2*lsum+1) * sqrt(pi / (gamma2+gamma)) 
                   / ( (gamma+gamma2)**(lsum+1) * 2**(lsum+2) ))
        cff = self._coef(lsum, lxx, lyy, lzz, gamma2) * self._coef(lsum, lx, ly, lz, gamma)
        c_inv_cmplx = (self._sphr_to_dec_cf(lsum, ml, lxx, lyy, lzz).conjugate()
                       * self._sphr_to_dec_cf(lsum, ml, lx, ly, lz))
        c_inv_real = c_inv_cmplx.real
        return cff * c_inv_real * rad_int

















                

