 The script is designed to process the Dirac program output file: extract all the necessary 
 information and put it into the file formatted to be recognized by the ElSA program.
 The Cartesian Gaussian basis set is read from the file provided with the Dirac program 
 (here the dyall.v4z file in the . directory, but in principle it may change depending of 
 the parameters of the Dirac calculations). 
 The coefficients of the 3D atomic wave functions are reduced to the coefficients of only 
 1D atomic radial wave functions under the assumption that the angular atomic parts are just 
 spherical harminics. At the intermediate step all atomic radial orbitals are also represented 
 on the radial grid in the direct space and written to the corresponding files for checking.

 All the *.py files should be place to some working directory. In this directory the results
 of the Dirac program runs should be placed to the different subdirectories named after the 
 corresponding atomic species (./atom). Thus each subdirectory which name consists of no more 
 than 2 symbols is considered as a Dirac program run subdirectory. In each such subdirectory the 
 file with the name "atom_atom.out" is sought. This file is the only file to be read for each 
 specific atomic species and it must contain the needed atomic wave function coefficients. 
 In order to make the Dirac program to produce a correct *.out file the following strings 
 should be put to the Dirac input file:

 **ANALYZE
 .PRIVEC
 *PRIVEC
 .PRICMP
 1 0

 Note: For the ElSA program we are interested in only large components of closed inactive 
 atomic orbitals. This is why all the small components as well as open shell and active 
 orbitals are neglected in this script as well.
 
 The relaticistic basis-set-containing file should be also provided and put to the working
 directory. These files are provided with the Dirac program. In this script it is named as
 ./dyall.v4z file. This name may be changed accordingly.

 An example of proper Dirac input and output files for the Al atom is placed to 
 the ./Al-example directory. In order to run the example just rename the directory to ./Al and run 
 the main.py file.

 In case of any questions please contact me under r.e.ovcharenko@gmail.com
