import numpy as np
from numpy import genfromtxt
import random as rnd

N1UP = genfromtxt(r"numerical results/CP Interpolation/N1_Up.csv", delimiter=',')
N1DOWN = genfromtxt(r"numerical results/CP Interpolation/N1_Down.csv", delimiter=',')

N2UP = genfromtxt(r"numerical results/CP Interpolation/N2_Up.csv", delimiter=',')
N2DOWN = genfromtxt(r"numerical results/CP Interpolation/N2_Down.csv", delimiter=',')

N0UP = genfromtxt(r"numerical results/CP Interpolation/N0_Up.csv", delimiter=',')
N0DOWN = genfromtxt(r"numerical results/CP Interpolation/N0_Down.csv", delimiter=',')

NAVEUP = np.zeros(N1UP.shape)
NAVEDOWN = np.zeros(N1UP.shape)

for i in range(101):
    for j in range(101):
        n = rnd.randint(1, 3)
        if n == 1:
            NAVEUP[i, j] = N1UP[i, j]
            NAVEDOWN[i, j] = N1DOWN[i, j]
        if n == 2:
            NAVEUP[i, j] = N0UP[i, j]
            NAVEDOWN[i, j] = N0DOWN[i, j]
        if n == 3:
            NAVEUP[i, j] = N2UP[i, j]
            NAVEDOWN[i, j] = N2DOWN[i, j]

np.savetxt("NAVEUP.csv", NAVEUP, delimiter=",")
np.savetxt("NAVEDOWN.csv", NAVEDOWN, delimiter=",")
