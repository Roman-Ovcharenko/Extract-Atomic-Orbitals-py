##################################################################################################
# This file contains some general functions and subroutines which are not attributed to 
# a specific class
##################################################################################################
from math import sqrt
from math import exp
from math import factorial as fact
from math import ceil
from math import floor

##################################################################################################
# Print elapsed time to the screen 
##################################################################################################
def print_time(time):
    hours = int(time // 3600)
    remain = time - 3600 * hours
    minutes = int(remain // 60)
    seconds = int(remain - 60 * minutes)
    print("Elapsed time: {} h {} min {} sec".format(hours, minutes, seconds))
    return None

##################################################################################################
# 1D numerical integration of multiplication of the f(R) and g(R) functions represented on 
# the R radial grid
##################################################################################################
def integrate(R, f, g):
    Nr = len(R)
    if len(f) != Nr or len(g) != Nr:
        raise Exception("f, g and dR have inconsistent size")
    res = 0.0
    for ir in range(0, Nr-2, 2):
        a = R[ir]
        b = R[ir+2]
        fa = f[ir] * g[ir]
        fb = f[ir+2] * g[ir+2]
        fc = f[ir+1] * g[ir+1]
        res = res + (fa + 4*fc + fb) * (b - a) / 6.0
    if (Nr-1)%2 == 1:
        a = R[Nr-2]
        b = R[Nr-1]
        fa = f[Nr-2] * g[Nr-2]
        fb = f[Nr-1] * g[Nr-1]
        res = res + (b - a) * (fa + fb) / 2.0
    return res

##################################################################################################
# Implementation of the double factorial 
##################################################################################################
def dfactorial(n):
    if n < -1:
        raise Exception("Double factorial of the number < -1 is not defined")
    if n == -1 or n == 0:
        return 1
    if n % 2 == 0:
        i1 = 2
    else:
        i1 = 1
    res = 1
    for i in range(i1, n+1, 2):
        res = res * i
    return res

##################################################################################################
# The Binomial coefficients
##################################################################################################
def binom(n, k):
    if n < 0:
        raise Exception("n cannot be negative")
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    k = min(k, n - k) # take advantage of symmetry
    c = 1
    for i in range(k):
        c = c * (n - i) / (i + 1)
    return int(c)

##################################################################################################
# The symbolic representation of the orbital symmetry
##################################################################################################
def get_symb_ll(ll):
    if ll == 0:
        return "s"
    elif ll == 1:
        return "p"
    elif ll == 2:
        return "d"
    elif ll == 3:
        return "f"
    elif ll == 4:
        return "g"
    elif ll == 5:
        return "h"
    elif ll == 6:
        return "i"
    elif ll >= 7:
        raise Exception("Symbol for ll >= 5 is not implemented")
    return None

##################################################################################################
# Convert the rational fraction into the decimal fraction
##################################################################################################
def convert_to_float(str_):
    try:
        float_ = float(str_)
    except ValueError:
        num, denom = str_.split("/")
        float_ = float(num) / float(denom)
    return float_

##################################################################################################
# Implementation of the Clebsch-Gordon coefficients
##################################################################################################
def clebsch(c,gama,a,alpha,b,beta):
    eps = float(1e-6)
    if abs(gama) > c or abs(alpha) > a or abs(beta) > b:
        return 0.0
    if abs(gama - alpha - beta) > eps or b < abs(c-a) or b > c+a: 
        return 0.0
    dabc = sqrt( fact(a+b-c) * fact(a-b+c) * fact(b+c-a) / fact(a+b+c+1) )
    cf = sqrt( fact(c+gama) * fact(c-gama) * (2.0*c+1.0) 
              / ( fact(a+alpha) * fact(a-alpha) * fact(b+beta) * fact(b-beta) ) )

    Zmin = ceil(max(0.0, alpha-a, b+gama-a))
    Zmax = floor(min( c+b+alpha, c-a+b, c+gama))
    summ = 0.0
    for z in range(Zmin, Zmax+1):
        summ = summ + ( (-1.0)**(b+beta+z) * fact(c+b+alpha-z) * fact(a-alpha+z) 
                       / ( fact(z) * fact(c-a+b-z) * fact(c+gama-z) * fact(a-b-gama+z) ))
    return dabc * cf * summ

##################################################################################################
# Parse a symbolic representation of the orbital symmetry
##################################################################################################
def parse_orb_symm(sym):
    if sym[0] == "s":
        lx = 0
        ly = 0
        lz = 0
        ltot = 0
        return lx, ly, lz, ltot
    elif sym[0] == "p":
        lx = sym.count("x")
        ly = sym.count("y")
        lz = sym.count("z")
        ltot = lx + ly + lz
        return lx, ly, lz, ltot
    elif sym[0] == "d":
        lx = sym.count("x")
        ly = sym.count("y")
        lz = sym.count("z")
        ltot = lx + ly + lz
        return lx, ly, lz, ltot
    elif sym[0] == "f":
        lx = sym.count("x")
        ly = sym.count("y")
        lz = sym.count("z")
        ltot = lx + ly + lz
        return lx, ly, lz, ltot
    elif sym[0] == "g":
        lx = int(sym[1])
        ly = int(sym[2])
        lz = int(sym[3])
        ltot = lx + ly + lz
        return lx, ly, lz, ltot
    elif sym[0] == "h":
        lx = int(sym[1])
        ly = int(sym[2])
        lz = int(sym[3])
        ltot = lx + ly + lz
        return lx, ly, lz, ltot
    elif sym[0] == "i":
        lx = int(sym[1])
        ly = int(sym[2])
        lz = int(sym[3])
        ltot = lx + ly + lz
        return lx, ly, lz, ltot
    else:
        raise Exception("Such orbital symmetry: {} was not implemented".format(sym[0]))
    return None

##################################################################################################
# Read the next nonempty line from the file
##################################################################################################
def get_line(f):
    while True:
        line = f.readline()
        if line.strip():
            return line
    raise Exception("End of file")
    return None

##################################################################################################
# Get the representation of the radial atomic orbital in the direct space
##################################################################################################
def reconstruct_rad_orbital(cf_rad, exps, Rc, ll):
    phi = [0.0 for ir in range(len(Rc))]
    for ir in range(len(Rc)):
        sum_ = 0.0
        for iexp in range(len(exps)):
            sum_ = sum_ + cf_rad[iexp] * exp(-exps[iexp] * Rc[ir]**2)
        phi[ir] = Rc[ir]**(ll+1) * sum_
    return phi


   
    



    




    








