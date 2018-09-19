#!/usr/local/bin/python3.7
import os
import time

import basis
import atom
import others

file_exp = "../dyall.v4z"

path = os.getcwd()
print("Original root directory: {}".format(path))
root_tot = os.listdir(".")
print("It contains:\n{}".format(root_tot))

root = list(root_tot)
for dir_ in root_tot:
    if len(dir_) > 2 or os.path.isfile(dir_):
        root.remove(dir_)
print("List of atoms is:\n{}".format(root))
print()

tot_time_begin = time.time()
for at_symb in root:
    time_begin = time.time()
    os.chdir(at_symb+"/")
    path = os.getcwd()
    print("{:=^80s}".format(" "+at_symb+" "))                           
    print("Path: {}".format(path))
    print()
    file_prop = at_symb + "_" + at_symb + ".out"
    file_coef = file_prop
    file_ElSA = at_symb + "-dyall-v4z.h"
    bas = basis.Basis(file_prop, file_exp, at_symb)
    bas.printout()
    at = atom.Atom(file_prop, file_coef, bas, at_symb)
    at.printout()
    at.print_ElSA_file(file_ElSA)
    os.chdir("..")
    time_end = time.time()
    others.print_time(time_end - time_begin)
print("{:=^80s}\n".format(" END "))                           
tot_time_end = time.time()
others.print_time(tot_time_end - tot_time_begin)




