#!/bin/bash
Z=13.
#############################
if [ "${2}" = "t" ]; then
	echo " TEST CALCULATION"
	     NUCMOD=1
             GAUNT="#.GAUNT"
             DOSSSS="#.DOSSSS"
             LVCORR=".LVCORR"
             URKBAL="#.URKBAL"
else
	     NUCMOD=2
             GAUNT=".GAUNT"
             DOSSSS=".DOSSSS"
             LVCORR="#.LVCORR"
             URKBAL="#.URKBAL"
fi
#############################
cat > ${1}.inp <<ENDF
**DIRAC
.TITLE
${1} atom
.WAVE FUNCTION
.ANALYZE
#################
**WAVE FUNCTION
.SCF
*SCF
.OPEN SHELL
1
3/2,6
.CLOSED SHELL
 4 6
.MAXITR
150
.EVCCNV
1.0D-8 1.0D-8
#################
**ANALYZE
.PRIVEC
*PRIVEC
.PRICMP
 1 0
#################
**GENERAL
.DIRECT
 1 1 1
.SPHTRA
 1 1
.PCMOUT
.LINDEP
1.0D-6 1.0D-8
.CVALUE
137.0359998
#################
**HAMILTONIAN
${GAUNT}
${DOSSSS}
${LVCORR}
${URKBAL}
.INTFLG
1 1 1
#################
**INTEGRALS
.NUCMOD
${NUCMOD}
*READIN
.UNCONTRACT
*TWOINT
.AOFOCK
*END OF INPUT
ENDF
#############################
cat > ${1}.mol <<ENDF
DIRAC
${1} ATOM
smallest basis set
C 1                A    
       ${Z}  1   
${1}     0.000000  0.000000  0.000000
LARGE BASIS dyall.v4z
FINISH
ENDF
#############################
/Users/campic/Dirac-programm/DIRAC-17.0-Source/build/pam --inp=${1}.inp --mol=${1}.mol --get "DFPCMO" --noarch
#############################
echo 
echo " DONE!"
