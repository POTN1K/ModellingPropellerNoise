import numpy as np

# Blade characteristics
rpm = 8000
omega = rpm * 2 * np.pi / 60
r = 0.15
v_tip = omega * r
v_free = 9
j = v_free / (rpm / 60 * 2 * r)


def chord_poly(x, r_=0.15):
    chord_ratio = 5.31329 - 87.1748*x + 608.234*x**2 - 2289.46*x**3 + 5136.01*x**4 - 7082.1*x**5 + 5891.59*x**6 - 2714.5*x**7 + 532.151*x**8
    return chord_ratio * r_


alpha = np.arctan(v_free / v_tip) * 180 / np.pi
#print(alpha)
