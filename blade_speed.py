import numpy as np
import pandas as pd


def chord_poly(x, r_=0.15):
    """Polynomial giving chord at position x. 'x' only from 0.33-1"""
    chord_ratio = 5.31329 - 87.1748 * x + 608.234 * x ** 2 - 2289.46 * x ** 3 + 5136.01 * x ** 4 - 7082.1 * x ** 5 + 5891.59 * x ** 6 - 2714.5 * x ** 7 + 532.151 * x ** 8
    return chord_ratio * r_


def twist_poly(x):
    l = -8971.3298 + 206093.7010 * x - 2.0708079035e6 * (x ** 2) + 1.210070e7 * (x ** 3) - 4.5793508e7 * (
            x ** 4) + 1.18060411e8 * (x ** 5) - 2.1196809e8 * (x ** 6) + 2.65500274e8 * (x ** 7) - 2.277350e8 * (
                x ** 8) + 1.2759e8 * (x ** 9) - 4.20825e7 * (x ** 10) + 6.1972e6 * (x ** 11)
    l = -8971.329845161235 + 206093.7010106019 * x - 2.0708079035287888e6 * x ** 2 + 1.2100709270279987e7 * x ** 3 - 4.5793508020059824e7 * x ** 4 + 1.1806041107129405e8 * x ** 5 - 2.1196809140788797e8 * x ** 6 + 2.6550027408285195e8 * x ** 7 - 2.2773507697170436e8 * x ** 8 + 1.2759430429713346e8 * x ** 9 - 4.208253179063628e7 * x ** 10 + 6.197207001092323e6 * x ** 11
    return l


# Blade characteristics
rpm = 8000
omega = rpm * 2 * np.pi / 60
v_free = 9*1.25
j = v_free / (rpm / 60 * 0.3)


def bladeSection(r):
    v_tip = omega * r
    v_tot = np.sqrt(v_tip ** 2 + v_free ** 2)
    alpha = np.arctan(v_free / v_tip) * 180 / np.pi
    x = r / 0.15
    t = twist_poly(x)
    phi = -alpha + t
    print(phi)
    return v_tot, x, phi


r = np.linspace(0.033, .15, 10)

v_list = []
x_list = []
phi_list = []
for element in r:
    v, x, p = bladeSection(element)
    v_list.append(v)
    x_list.append(x)
    phi_list.append(p)

df = pd.DataFrame({"r": r, "x": x_list, "v": v_list, "phi": phi_list})
df.to_csv(r"measurements/loads/blade25.csv", index=None)
