from others import get_line
from others import parse_orb_symm

class Basis(object):

    _bohr_AA = float(0.529177249)

    def __init__(self, file_prop, file_exp, at_symb):
        self.at_symb = at_symb
        sizes_L, sizes_symb_L = self._get_size("large", file_prop)
        sizes_S, sizes_symb_S = self._get_size("small", file_prop)
        self.exps = self._get_exp(file_exp)
        (self.size,  self.symm, self.orb_symm, self.exp, 
         self.lx, self.ly, self.lz, self.ltot) = self._gen_basis(file_prop, sizes_L, 
                                                             sizes_symb_L, sizes_S)

    def printout(self):
        print("BASIS OUTPUT")
        print("Basis set size (size): {}".format(self.size))
        print("The entire basis set has been written to the basis.dat file")
        with open("basis.dat", "w") as f:
            for ibas in range(self.size):
                f.write("{}   {}   {}   {}   {}   {}   {}   {}\n".
                        format(ibas, self.symm[ibas], self.orb_symm[ibas], self.lx[ibas], 
                               self.ly[ibas], self.lz[ibas], self.ltot[ibas], self.exp[ibas]))
        print("{}\n".format("-"*80))
        return None

    def _get_size(self, str_, file_):
        with open(file_) as fp:
            for line in fp:
                if (str_ == "large") and ("Number of large orbitals in each symmetry:" in line):
                    _, _, _, _, _, _, _, Ag, B3u, B2u, B1g, B1u, B2g, B3g, Au =  line.split()
                    break
                if (str_ == "small") and ("Number of small orbitals in each symmetry:" in line):
                    _, _, _, _, _, _, _, Ag, B3u, B2u, B1g, B1u, B2g, B3g, Au =  line.split()
                    break
            return ([int(Ag), int(B3u), int(B2u), int(B1g), int(B1u), int(B2g), 
                     int(B3g), int(Au)], 
                    ["Ag", "B3u", "B2u", "B1g", "B1u", "B2g", "B3g", "Au"])

    def _get_exp(self, file_):
        with open(file_) as f:
            for _ in range(29):
                line = f.readline()
            at_char = line[1:].strip()
            while self.at_symb != at_char:
                line = f.readline()
                at_char = line[1:].strip()
            f.readline()
            exp_s = self._get_exp_sym(f, " s ")
            exp_p = self._get_exp_sym(f, " p ")
            exp_d = self._get_exp_sym(f, " d ")
            exp_f = self._get_exp_sym(f, " f ")
            exp_g = self._get_exp_sym(f, " g ")
            exp_h = self._get_exp_sym(f, " h ")
            exp_i = self._get_exp_sym(f, " i ")

        return [exp_s, exp_p, exp_d, exp_f, exp_g, exp_h, exp_i]

    def _gen_basis(self,file_, sizes_L, sizes_symb_L, sizes_S):
        size = 0
        exp_gen = []
        sym_gen = []
        orb_sym_gen = []
        lx_gen = []
        ly_gen = []
        lz_gen = []
        ltot_gen = []
        qu_sym_gen = []
        qu_orb_sym_gen = []
        qu_lx_gen = []
        qu_ly_gen = []
        qu_lz_gen = []
        qu_ltot_gen = []
        with open(file_) as f:
            for line in f:
                if "Large component functions" in line:
                    break
            line = get_line(f)
            for isym in range(len(sizes_L)):
                if sizes_symb_L[isym] in line:
                    ltot_prev = -10
                    queue = 0
                    while True:
                        line = get_line(f)
                        if ("Symmetry" in line) or ("Small component functions" in line):
                            size, queue = self._update_counters(queue, size, sym_gen, 
                                                                qu_sym_gen, orb_sym_gen, 
                                                                qu_orb_sym_gen, lx_gen, 
                                                                qu_lx_gen, ly_gen, qu_ly_gen,
                                                                lz_gen, qu_lz_gen, ltot_gen,
                                                                qu_ltot_gen, exp_gen, 
                                                                self.exps[ltot_prev])
                            break
                        lst = line.split()
                        num_loc = int(lst[0])
                        sym = lst[-1]
                        lx, ly, lz, ltot = parse_orb_symm(sym)
                        if num_loc != len(self.exps[ltot]):
                            raise Exception("The number of exponents from " 
                                            "the basis-set-containing-file: {} \n"
                                            "does not correspond to the numbers in "
                                            "the dirac output file: {}"
                                            .format(len(self.exps[ltot]), num_loc))
                        if (ltot == ltot_prev) or (queue == 0):
                            queue = self._append_to_counters(qu_sym_gen, qu_orb_sym_gen,
                                                             qu_lx_gen, qu_ly_gen, qu_lz_gen,
                                                             qu_ltot_gen, sizes_symb_L[isym], 
                                                             sym, lx, ly, lz, ltot, queue)
                        else:
                            size, queue = self._update_counters(queue, size, sym_gen, 
                                                                qu_sym_gen, orb_sym_gen, 
                                                                qu_orb_sym_gen, lx_gen, 
                                                                qu_lx_gen, ly_gen, qu_ly_gen, 
                                                                lz_gen, qu_lz_gen, ltot_gen,
                                                                qu_ltot_gen, exp_gen, 
                                                                self.exps[ltot_prev])
                            queue = self._append_to_counters(qu_sym_gen, qu_orb_sym_gen,
                                                             qu_lx_gen, qu_ly_gen, qu_lz_gen,
                                                             qu_ltot_gen, sizes_symb_L[isym], 
                                                             sym, lx, ly, lz, ltot, queue)
                        ltot_prev = ltot
                else:
                    raise Exception("Symmetry symbol: {} is not found".format(sizes_symb_L[isym]))
                size = size + sizes_S[isym]
                for i in range(sizes_S[isym]):
                    exp_gen.append(-1.0)
                    sym_gen.append("S")
                    orb_sym_gen.append("S")
                    lx_gen.append(-1)
                    ly_gen.append(-1)
                    lz_gen.append(-1)
                    ltot_gen.append(-1)

        return size,  sym_gen, orb_sym_gen, exp_gen, lx_gen, ly_gen, lz_gen, ltot_gen

    def _get_exp_sym(self, f, sym):
        exp = []
        line = f.readline()
        if sym in line:
            line = f.readline()
            line_list = line.split()
            for i in range(int(line_list[0])):
                line = f.readline()
                exp.append(float(line) / self._bohr_AA**2)
        return exp

    def _update_counters(self, queue, num, sym_gen, qu_sym_gen, orb_sym_gen, qu_orb_sym_gen,
                         lx_gen, qu_lx_gen, ly_gen, qu_ly_gen, lz_gen, qu_lz_gen, ltot_gen,
                         qu_ltot_gen, exp_gen, exp):
        for iexp in range(len(exp)):
            num = num + queue
            sym_gen.extend(qu_sym_gen)
            orb_sym_gen.extend(qu_orb_sym_gen)
            lx_gen.extend(qu_lx_gen)
            ly_gen.extend(qu_ly_gen)
            lz_gen.extend(qu_lz_gen)
            ltot_gen.extend(qu_ltot_gen)
            for iqu in range(queue):
                exp_gen.append(float(exp[iexp]))
        queue = 0
        qu_sym_gen.clear()
        qu_orb_sym_gen.clear()
        qu_lx_gen.clear()
        qu_ly_gen.clear()
        qu_lz_gen.clear()
        qu_ltot_gen.clear()
        return num, queue

    def _append_to_counters(self, qu_sym_gen, qu_orb_sym_gen, qu_lx_gen, qu_ly_gen, qu_lz_gen,
                            qu_ltot_gen, sizes, sym, lx, ly, lz, ltot, queue):
        queue = queue + 1
        qu_sym_gen.append(sizes)
        qu_orb_sym_gen.append(sym)
        qu_lx_gen.append(lx)
        qu_ly_gen.append(ly)
        qu_lz_gen.append(lz)
        qu_ltot_gen.append(ltot)
        return queue
 




