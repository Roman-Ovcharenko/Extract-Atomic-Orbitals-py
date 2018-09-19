import orbital
import red_orbital
from others import convert_to_float
from others import dfactorial as dfact

class Atom(object):

    def __init__(self, file_prop, file_coef, bas, at_symb):
        self.file_prop = file_prop
        self.file_coef = file_coef
        self.at_symb = at_symb
        self.bas = bas
        (energies, jjs, lls, mjs, nns, parities, 
         self.num_orb_ger, self.num_orb_ung, self.num_orb, cfs_La, cfs_Lb) = self._get_param()
        self.orbitals = []
        for iorb in range(self.num_orb):
            self.orbitals.append(orbital.Orbital(self.file_coef, energies[iorb], parities[iorb], 
                                                 jjs[iorb], mjs[iorb], lls[iorb], self.bas, iorb, 
                                                 self.num_orb_ger, self.at_symb, cfs_La[iorb],
                                                 cfs_Lb[iorb], nns[iorb]))
        orb_avg = self._average_over_mj()
        self.num_red_orb = len(orb_avg)
        self.red_orbitals = []
        for irorb in range(self.num_red_orb):
            self.red_orbitals.append(red_orbital.RedOrbital(self.orbitals, irorb, 
                                                            orb_avg[irorb], self.at_symb, bas))
        return None

    def printout(self):
        print("\nATOM OUTPUT")
        print("Total number of orbitals: {}".format(self.num_orb))
        print("Number of gerade orbitals: {}".format(self.num_orb_ger))
        print("Number of ungerade orbitals: {}\n".format(self.num_orb_ung))
        print("iorb       Symbol         nn    ll    jj      mj       E / hartree")
        for iorb in range(self.num_orb):
            print("{:>3d}    {:<14s}    {:>3d}    {:>2d}    {:>3.1f}    {:>+4.1f}    {:>11.3f}"
                  .format(iorb, self.orbitals[iorb].orb_sym, self.orbitals[iorb].nn, 
                          self.orbitals[iorb].ll, self.orbitals[iorb].jj, self.orbitals[iorb].mj,
                          self.orbitals[iorb].energy))
        print()
        print("Total number of the orbitals reduced by averaging over mj: {}\n"
              .format(self.num_red_orb))
        print("iorb    Symbol         nn    ll    jj       E / hartree      Averaged over (mj)")
        for irorb in range(self.num_red_orb):
            mj_list = ([self.orbitals[self.red_orbitals[irorb].avg_list[jorb]].mj 
                        for jorb in range(self.red_orbitals[irorb].num_deg_orb)])
            print("{:>3d}    {:<14s} {:>3d}    {:>2d}    {:>3.1f}    {:>11.3f}        {}"
                  .format(irorb, self.red_orbitals[irorb].orb_sym, self.red_orbitals[irorb].nn, 
                          self.red_orbitals[irorb].ll, self.red_orbitals[irorb].jj,
                          self.red_orbitals[irorb].energy,
                          mj_list))
        print("{}\n".format("-"*80))
        return None

    def _get_ll(self, jj, parities):
        ll = []
        for iorb in range(self.num_orb):
            l_minus = int(jj[iorb] - 0.5)
            l_plus = int(jj[iorb] + 0.5)
            if parities[iorb] == "gerade" and l_minus % 2 == 0:
                ll.append(l_minus)
            elif parities[iorb] == "gerade" and l_plus % 2 == 0:
                ll.append(l_plus)
            elif parities[iorb] == "ungerade" and l_minus % 2 != 0:
                ll.append(l_minus)
            elif parities[iorb] == "ungerade" and l_plus % 2 != 0:
                ll.append(l_plus)
        return ll

    def print_ElSA_file(self, file_):
        ll_max = self._get_ll_max()
        f = open(file_, "w+")
        f.write("!========== Element specific information =================\n")
        f.write('   sym_at = "{:>3s}"\n'.format(self.at_symb))
        f.write('   Norb = {}\n'.format(self.num_red_orb))

        str_ = "   Nbas(1:Nsym) = (/ {}".format(len(self.bas.exps[0]))
        for il in range(1, ll_max+1):
            str_ = str_ + ", {}".format(len(self.bas.exps[il]))
        str_ = str_ + " /)\n"
        f.write(str_)

        f.write('#include "../getRelOrb.h"\n')
        f.write("!==================== Atomic orbital's parameters ============================\n")
        str_orb_sym = "   orbstr(1:Norb) = (/ &\n"
        str_en = "   orbE(1:Norb) = (/ &\n"
        str_nn = "   nn(1:Norb) = (/ &\n"
        str_ll = "   ll(1:Norb) = (/ &\n"
        str_jj = "   jj(1:Norb) = (/ &\n"
        for iorb in range(self.num_red_orb-1):
            str_orb_sym = str_orb_sym + '    "{}", &\n'.format(self.red_orbitals[iorb].orb_sym)
            str_en = str_en + ' {:>18.10f}, &\n'.format(self.red_orbitals[iorb].energy)
            str_nn = str_nn + '      {}, &\n'.format(self.red_orbitals[iorb].nn)
            str_ll = str_ll + '      {}, &\n'.format(self.red_orbitals[iorb].ll)
            str_jj = str_jj + ' {:6.1f}, &\n'.format(self.red_orbitals[iorb].jj)
        str_orb_sym = (str_orb_sym 
                       + '    "{}" /)\n\n'.format(self.red_orbitals[self.num_red_orb-1].orb_sym))
        str_en = (str_en + " {:>18.10f} /) ! hartree\n\n"
                  .format(self.red_orbitals[self.num_red_orb-1].energy))
        str_nn = str_nn + "      {} /)\n\n".format(self.red_orbitals[self.num_red_orb-1].nn)
        str_ll = str_ll + "      {} /)\n\n".format(self.red_orbitals[self.num_red_orb-1].ll)
        str_jj = str_jj + " {:6.1f} /)\n\n".format(self.red_orbitals[self.num_red_orb-1].jj)

        f.write(str_orb_sym)
        f.write(str_en)
        f.write(str_nn)
        f.write(str_ll)
        f.write(str_jj)

        f.write("!====================== Basis exponents ============================\n")

        for il in range(ll_max+1):
            f.write("   isym = {}\n".format(il+1))
            f.write("   alpha(1:Nbas(isym),isym) = (/ &\n")
            for iexp in range(len(self.bas.exps[il])-1):
                f.write("{:>26.16e}, &\n".format(self.bas.exps[il][iexp]))
            f.write("{:>26.16e} /)\n\n".format(self.bas.exps[il][len(self.bas.exps[il])-1]))

        f.write("!====================== Atomic wave functions  ============================\n")
        for iorb in range(self.num_red_orb):
            f.write("   iorb = {:d}\n".format(iorb+1))
            f.write("   cf( 1: , iorb ) = (/ &\n")
            for iexp in range(len(self.bas.exps[self.red_orbitals[iorb].ll])-1):
                f.write("{:>26.16e}, &\n".format(self.red_orbitals[iorb].cf[iexp]))
            f.write("{:>26.16e} /)\n\n".format(self.red_orbitals[iorb].
                                         cf[len(self.bas.exps[self.red_orbitals[iorb].ll])-1]))
        f.write("!====================== End of element specific block "
                "============================\n")
        f.close()
        return None

    def _get_param(self):
        energies_ger = []
        energies_ung = []
        with open(self.file_prop) as f:
            line = f.readline()
            while line:
                line = f.readline()
                if "Occupation in fermion symmetry E1g" in line:
                    num_orb_ger, mjs_ger, parities_ger = self._get_num_orb(f, "gerade")
                elif "Occupation in fermion symmetry E1u" in line:
                    num_orb_ung, mjs_ung, parities_ung = self._get_num_orb(f, "ungerade")

                if "Fermion ircop E1g" in line: 
                    energies_ger = self._get_energies(f, num_orb_ger)
                if "Fermion ircop E1u" in line: 
                    energies_ung = self._get_energies(f, num_orb_ung)

        cfs_La_ger = []
        cfs_Lb_ger = []
        lls_ger = []
        for iorb in range(num_orb_ger):
            La_ger, Lb_ger = self._get_coef(parities_ger[iorb], iorb, energies_ger[iorb])
            cfs_La_ger.append(La_ger)
            cfs_Lb_ger.append(Lb_ger)
            lls_ger.append(self._get_ll(parities_ger[iorb], cfs_La_ger[iorb], cfs_Lb_ger[iorb]))

        cfs_La_ung = []
        cfs_Lb_ung = []
        lls_ung = []
        for iorb in range(num_orb_ung):
            La_ung, Lb_ung = self._get_coef(parities_ung[iorb], iorb, energies_ung[iorb])
            cfs_La_ung.append(La_ung)
            cfs_Lb_ung.append(Lb_ung)
            lls_ung.append(self._get_ll(parities_ung[iorb], cfs_La_ung[iorb], cfs_Lb_ung[iorb]))

        jjs_ger = self._get_jj(mjs_ger, energies_ger, lls_ger)
        jjs_ung = self._get_jj(mjs_ung, energies_ung, lls_ung)

        nns_ger = self._get_nn(jjs_ger, lls_ger, mjs_ger, energies_ger)
        nns_ung = self._get_nn(jjs_ung, lls_ung, mjs_ung, energies_ung)

        num_orb = num_orb_ger + num_orb_ung
        jjs = jjs_ger + jjs_ung
        lls = lls_ger + lls_ung
        mjs = mjs_ger + mjs_ung
        nns = nns_ger + nns_ung
        self._assign_sign(mjs)
        parities = parities_ger + parities_ung
        energies = energies_ger + energies_ung
        cfs_La = cfs_La_ger + cfs_La_ung
        cfs_Lb = cfs_Lb_ger + cfs_Lb_ung

        return (energies, jjs, lls, mjs, nns, parities, num_orb_ger, 
                num_orb_ung, num_orb, cfs_La, cfs_Lb)

    def _get_nn(self, jjs, lls, mjs, energies):
        if len(mjs) == 0:
            return []
        _eps = 1.0e-6
        nns = []
        del_nn = [0] * len(mjs)
        for iorb in range(len(mjs)):
            for jorb in range(len(mjs)):
                if (abs(jjs[iorb] - jjs[jorb]) < _eps 
                    and lls[iorb] == lls[jorb] 
                    and abs(mjs[iorb] - mjs[jorb]) < _eps
                    and energies[jorb] < energies[iorb]):
                    del_nn[iorb] = del_nn[iorb] + 1
            nns.append(lls[iorb] + 1 + del_nn[iorb])
        return nns

    def _get_ll(self, parity, cf_La, cf_Lb):
        if parity == "gerade":
            ll_list = [0, 2, 4]
        else:
            ll_list = [1, 3, 5]
        contrib_max = 0.0
        ll_cur = -1
        for il, ll in enumerate(ll_list):
            contrib = 0.0
            for beta in [-0.5, 0.5]:
                contrib = contrib + self._get_contrib(ll, beta, cf_La, cf_Lb)
            if contrib > contrib_max:
                contrib_max = contrib
                ll_cur = ll
        if ll_cur < 0:
            raise Exception("ll_cur was not assigned")
        return ll_cur

    def _get_contrib(self, ltot, beta, cf_La, cf_Lb): 
        eps = 1.0e-10
        if abs(beta - 0.5) < eps:
            cf = cf_La
        else:
            cf = cf_Lb
        sum_ = 0.0
        for ibas in range(self.bas.size):
            if self.bas.ltot[ibas] != ltot:
                continue
            dens_mx = abs(cf[ibas])**2
            overlap = (dfact(2*self.bas.lx[ibas]-1) * dfact(2*self.bas.ly[ibas]-1) 
                       * dfact(2*self.bas.lz[ibas]-1))
            sum_ = sum_ + dens_mx * overlap
        return sum_

    def _get_coef(self, parity, iorb, energy):
        epsilon = float(2.0e-11)
        La = [complex(0.0, 0.0) for ibas in range(self.bas.size)]
        Lb = [complex(0.0, 0.0) for ibas in range(self.bas.size)]
        ival = -1
        nval = iorb 

        with open(self.file_coef) as f:
            while True:
                line = f.readline()
                if (parity == "gerade") and ("Fermion ircop E1g" in line):
                    break
                if (parity == "ungerade") and ("Fermion ircop E1u" in line):
                    break

            while True:
                line = f.readline()
                if "Electronic eigenvalue no." in line:
                    ival = ival + 1
                if ival == nval:
                    break

            lst = line.split()
            if abs(float(lst[-1]) - energy) > epsilon:
                raise Exception("Energy: {} is different from: {} for iorb: {}"
                                .format(float(lst[-1]),  energy, iorb))

            line = f.readline()
            while True:
                line = f.readline()
                lst = line.split()
                try:
                    ibas = int(lst[0]) - 1
                except:
                    break

                rLa = float(lst[-4])
                iLa = float(lst[-3])
                rLb = float(lst[-2])
                iLb = float(lst[-1])

                if (self.bas.symm[ibas] in line) and (self.bas.orb_symm[ibas] in line):
                    La[ibas] = complex(rLa, iLa)
                    Lb[ibas] = complex(rLb, iLb)
                else:
                    raise Exception("The order of basis functions does not coincide with "
                                    "that of wave functions")

        return La, Lb

    def _assign_sign(self, mjs):
        for i in range(len(mjs)):
            mjs[i] = mjs[i] * (-1.0)**int(mjs[i]-0.5)
        return None

    def _get_num_orb(self, f, parity):
        line_sum = ""
        line = f.readline()
        if "Inactive orbitals" in line:
            while True:
                line = f.readline()
                if "orbitals" not in line: 
                    line_sum = line_sum + line
                else:
                    break
            mjs = [convert_to_float(i) for i in line_sum.split()]
            num_orb = len(mjs)
            parities = [parity for i in range(num_orb)]
            return num_orb, mjs, parities
        else:
            return 0, [], []

    def _get_energies(self, f, num_orb):
        energies = []
        line = f.readline()
        num = 0
        while num < num_orb:
            line = f.readline()
            if "Electronic eigenvalue no." in line:
                lst = line.split()
                energies.append(float(lst[-1]))
                num = num + 1
        return energies

    def _get_jj(self, mjs, energies, lls):
        if len(mjs) == 0:
            return []
        _eps = 1.0e-6
        jjs = [0.0 for i in mjs]
        if_taken = [False for i in mjs]
        jj_max = max(mjs)
        if (jj_max - 0.5 < _eps):
            jj_cur = 0.5
            for i in range(len(mjs)):
                if self._if_compatible_ll(jj_cur, lls[i]):
                    self._change_jj(jjs, if_taken, jj_cur, i)
        else:
            jj_cur = jj_max
            while True:
                if jj_cur <= 0.0:
                    break
                eng_prev, index_prev = self._get_energy(jj_cur, mjs, if_taken, 
                                                                energies, lls)
                self._change_jj(jjs, if_taken, jj_cur, index_prev)
                mj_cur = jj_cur - 1.0
                while True:
                    if mj_cur < 0.0:
                        break
                    eng_cur, index_cur, if_found = self._get_closest_energy(jj_cur, mj_cur,
                                                                            if_taken, energies,
                                                                            mjs, eng_prev, lls)
                    if not if_found:
                        mj_cur = mj_cur - 1.0
                        continue
                    self._change_jj(jjs, if_taken, jj_cur, index_cur)
                    eng_prev = eng_cur
                    mj_cur = mj_cur - 1.0
                jj_cur = self._decrease_jj_cur(jj_cur, mjs, if_taken, lls)

        for i in range(len(if_taken)):
            if not if_taken[i]:
                jj_cur = lls[i] - 0.5
                self._change_jj(jjs, if_taken, jj_cur, i)
        return jjs

    def _get_energy(self, jj_cur, mjs, if_taken, energies, lls) :
        _eps = 1.0e-6
        for i in range(len(mjs)-1, -1, -1):
            if ((abs(abs(mjs[i]) - jj_cur) < _eps) 
                and (not if_taken[i])
                and self._if_compatible_ll(jj_cur, lls[i])):
                eng = energies[i]
                index = i
                return eng, index
        raise Exception("All energies for given jj_cur are already taken")
        return None

    def _if_compatible_ll(self, jj_cur, ll):
        _eps = 1.0e-6
        if ((abs(jj_cur - (ll + 0.5)) < _eps) 
            or (abs(jj_cur - (ll - 0.5)) < _eps)):
            return True
        else:
            return False

    def _get_closest_energy(self, jj_cur, mj_cur, if_taken, energies, mjs, eng_prev, lls):
        _eps = 1.0e-6
        eps_eng = 1000.0
        index = -1
        eng_index = 0.0
        for i in range(len(mjs)-1, -1, -1):
            eng_tst = energies[i]
            if ((abs(mjs[i] - mj_cur) < _eps) 
                and (not if_taken[i]) 
                and (abs(eng_tst - eng_prev) < eps_eng)
                and (self._if_compatible_ll(jj_cur, lls[i]))):
                eps_eng = abs(eng_tst - eng_prev)
                eng_index = eng_tst
                index = i
        if index < 0:
            if_found = False
        else:
            if_found = True
        return eng_index, index, if_found

    def _change_jj(self, jjs, if_taken, jj_cur, index):
        jjs[index] = jj_cur
        if_taken[index] = True
        return None

    def _decrease_jj_cur(self, jj_cur, mjs, if_taken, lls):
        _eps = 1.0e-6
        for i in range(len(mjs)-1, -1, -1):
            if ((abs(abs(mjs[i]) - jj_cur) < _eps) 
                and (not if_taken[i])
                and self._if_compatible_ll(jj_cur, lls[i])):
                return jj_cur
        return jj_cur - 1.0

    def _average_over_mj(self):
        orb_avg = []
        if_taken = [False for iorb in range(self.num_orb)]
        for iorb in range(self.num_orb):
            if if_taken[iorb]:
                continue
            orb_cur_avg = []
            orb_cur = self.orbitals[iorb]
            if_taken[iorb] = True
            orb_cur_avg.append(iorb)
            for jorb in range(iorb, self.num_orb):
                if (not if_taken[jorb] 
                    and (orb_cur.nn == self.orbitals[jorb].nn) 
                    and (orb_cur.ll == self.orbitals[jorb].ll)
                    and (orb_cur.jj == self.orbitals[jorb].jj)
                    and (orb_cur.mj != self.orbitals[jorb].mj)):
                    orb_cur_avg.append(jorb)
                    if_taken[jorb] = True
            orb_avg.append(orb_cur_avg)
        return orb_avg

    def _get_ll_max(self):
        ll_max = -1
        for iorb in range(self.num_red_orb):
            ll = self.red_orbitals[iorb].ll
            if ll > ll_max:
                ll_max = ll
        if ll_max < 0:
            raise Exception("Wrong ll_max value")
        return ll_max














