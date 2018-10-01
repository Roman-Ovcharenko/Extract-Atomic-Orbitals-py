#!/usr/local/bin/python3.7
#
##################################################################################################
# The script is designed to process the Dirac program output file: extract all the necessary 
# information and put it into the file formatted to be recognized by the ElSA program.
# The Cartesian Gaussian basis set is read from the file provided with the Dirac program 
# (here the dyall.v4z file in the . directory, but in principle it may change depending of 
# the parameters of the Dirac calculations). 
# The coefficients of the 3D atomic wave functions are reduced to the coefficients of only 
# 1D atomic radial wave functions under the assumption that the angular atomic parts are just 
# spherical harminics. At the intermediate step all atomic radial orbitals are also represented 
# on the radial grid in the direct space and written to the corresponding files for checking.
#
# All the *.py files should be place to some working directory. In this directory the results
# of the Dirac program runs should be placed to the different subdirectories named after the 
# corresponding atomic species (./atom). Thus each subdirectory which name consists of no more 
# than 2 symbols is considered as a Dirac program run subdirectory. In each such subdirectory the 
# file with the name "atom_atom.out" is sought. This file is the only file to be read for each 
# specific atomic species and it must contain the needed atomic wave function coefficients. 
# In order to make the Dirac program to produce a correct *.out file the following strings 
# should be put to the Dirac input file:
#
# **ANALYZE
# .PRIVEC
# *PRIVEC
# .PRICMP
# 1 0
#
# Note: For the ElSA program we are interested in only large components of closed inactive 
# atomic orbitals. This is why all the small components as well as open shell and active 
# orbitals are neglected in this script as well.
# 
# The relaticistic basis-set-containing file should be also provided and put to the working
# directory. These files are provided with the Dirac program. In this script it is named as
# ./dyall.v4z file. This name may be changed accordingly.
#
# An example of proper Dirac input and output files for the Al atom is placed to 
# the ./Al-example directory. In order to run the example just rename the directory to ./Al and run 
# the main.py file.
#
# In case of any questions please contact me under r.e.ovcharenko@gmail.com
##################################################################################################
import os
import time

import basis
import atom
import others

# Define the file containing relativistic basis set used in the Dirac program
file_exp = "../dyall.v4z"

# Define the current directory
path = os.getcwd()
print("Original root directory: {}".format(path))
root_tot = os.listdir(".")
print("It contains:\n{}".format(root_tot))

# Define a list of atoms in the current directory
root = list(root_tot)
for dir_ in root_tot:
    if len(dir_) > 2 or os.path.isfile(dir_):
        root.remove(dir_)
print("List of atoms is:\n{}".format(root))
print()

# Run a loop over different atomic species
wrong_norm_orbitals = []
corr_norm_orbitals_symb = [] 
corr_norm_orbitals_val = []
lls_pre = [] 
dev_ll_contrib = []
dev_ll_contrib_symb = []
tot_time_begin = time.time()
for at_symb in root:
    time_begin = time.time()
    os.chdir(at_symb+"/")
    path = os.getcwd()
    print("{:=^80s}".format(" "+at_symb+" "))                           
    print("Path: {}".format(path))
    print()
    file_prop = at_symb + "_" + at_symb + ".out"
    file_ElSA = at_symb + "-dyall-v4z.h"
    bas = basis.Basis(file_prop, file_exp, at_symb)
    bas.printout()
    at = atom.Atom(file_prop, bas, at_symb, wrong_norm_orbitals, corr_norm_orbitals_symb,
                   corr_norm_orbitals_val, lls_pre, dev_ll_contrib, dev_ll_contrib_symb)
    at.printout()
    at.print_ElSA_file(file_ElSA)
    os.chdir("..")
    time_end = time.time()
    others.print_time(time_end - time_begin)
print("{:=^80s}\n".format(" END "))   

# Report errors and warnings occured during execution
print("The number of orbitals with wrong norm: {}:".format(len(wrong_norm_orbitals)))
for orb in wrong_norm_orbitals:
    print(orb)
print()
max_el = max(corr_norm_orbitals_val)
index_max = corr_norm_orbitals_val.index(max_el)
orb_max_el = corr_norm_orbitals_symb[index_max]
print("The maximum relative deviation of the norm was observed for the {} orbitals and it is: {} %"
        .format(orb_max_el, max_el))
print()
max_dev_ll_val = max(dev_ll_contrib)
index_max = dev_ll_contrib.index(max_dev_ll_val)
orb_max_dev_ll = dev_ll_contrib_symb[index_max]
ll_pre_val = lls_pre[index_max]
print("Maximum relative deviation in the assignment of the l quantum number was found "
      "for the {} orbital, equals to {:5.3f} % and corresponds to the l = {} orbital quantum number"
      .format(orb_max_dev_ll, max_dev_ll_val, ll_pre_val))
print()
tot_time_end = time.time()
others.print_time(tot_time_end - tot_time_begin)




