from urllib.request import urlopen
import ast
import os
import json
import numpy as np
import math as m
import csv
from math import sqrt,pi,cos,sin

def split_elems(formula):
    """Takes a list and splits strings inside it into title()'d pieces.
    Replaces the former string with the split stings:
        ['H',2,'SeO',4] -> ['H',2,'Se','O',4]"""
    for idx, name in enumerate(formula[:]):
        if isinstance(name,str):
            if sum(c.isupper() for c in name)>1:
                tmp = []
                for c in name:
                    if c.isupper():
                        tmp.append([c])
                    else:
                        tmp[-1].append(c)
                formula.pop(idx)
                for t in tmp[::-1]:
                    formula.insert(idx,"".join(t))
    return(formula)

def F2C(U):
    # Accepts: an arrays of the six unit cell parameters a,b,c,alpha,beta,gamma
    # Return: Transformation matrix from Fractional to Cartesian Coordinates
    rad = pi/180.
    a   = U[0]
    b   = U[1]
    c   = U[2]
    alp = rad*U[3]
    bet = rad*U[4]
    gam = rad*U[5]
    V   = a*b*c*sqrt(1-cos(gam)**2-cos(alp)**2-cos(bet)**2+2*cos(alp)*cos(bet)*cos(gam))
    M   = np.array([[a, b*cos(gam), c*cos(bet)],
        [0., b*sin(gam), c*(cos(alp)-cos(bet)*cos(gam))/(sin(gam))],
        [0., 0., V/(a*b*sin(gam))]])

    eps = 1E-10
    M[np.abs(M)<eps] = 0
    return(M)

def C2F(U):
    # Accepts: an arrays of the six unit cell parameters a,b,c,alpha,beta,gamma
    # Return: Transformation matrix from Cartesian to Fractional Coordinates
    rad = pi/180.
    a   = U[0]
    b   = U[1]
    c   = U[2]
    alp = rad*U[3]
    bet = rad*U[4]
    gam = rad*U[5]
    V   = a*b*c*sqrt(1-cos(gam)**2-cos(alp)**2-cos(bet)**2+2*cos(alp)*cos(bet)*cos(gam))
    M1  = np.array([[1/a,-cos(gam)/(a*sin(gam)),b*c*(cos(alp)*cos(gam)-cos(bet))/(V*sin(gam))],
          [0.,1/(b*sin(gam)),a*c*(cos(bet)*cos(gam)-cos(alp))/(V*sin(gam))],
          [0.,0.,a*b*sin(gam)/V]])

    eps = 1E-10
    M1[np.abs(M1)<eps] = 0
    return(M1)

def unique(list):
    x = []
    for a in list:
        if a not in x:
            x.append(a)
    x.sort()
    return(x)
