import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def txt_to_csv(file, name):
    read_file = pd.read_csv(file)
    read_file.to_csv("measurements/loads/" + name + ".csv", index=None)


def getForces():
    file = pd.read_csv("measurements/loads/no_grid.csv", header=None, usecols=[2, 5])
    # Bias that needs to be substracted
    bias = pd.read_csv("measurements/loads/bias.csv", header=None, usecols=[2, 5])
    print(bias.size)
    bias = bias.head(round(file.size / 2))
    print(bias.size)
    df = file.sub(bias)
    n = round(file.size / 2)

    x = pd.DataFrame(np.linspace(0, 20, n))

    df = pd.concat([df, x], axis=1, ignore_index=True)
    df.columns = ["Thrust", "Torque", "Time"]
    df.to_csv("measurements/loads/no_gridUPDATED.csv", index=None)


df = pd.read_csv("measurements/loads/no_gridUPDATED.csv")

fig, (ax1, ax2, ax3) = plt.subplots(3)
ax1.scatter(df["Time"], df["Thrust"])
ax2.scatter(df["Time"], df["Torque"])
ax3.scatter(df["Torque"], df["Thrust"])
plt.show()

Thrust_rms = np.sqrt(sum(df["Thrust"] ** 2) / df["Thrust"].size)
print(Thrust_rms)

Torque_rms = np.sqrt(sum(df["Torque"] ** 2) / df["Torque"].size)
print(Torque_rms)


def cx(cl, cd, phi):
    """Input phi in degrees, everything as numpy array"""
    phi = phi * np.pi / 180
    return cl * np.sin(phi) - cd * np.cos(phi)


def cy(cl, cd, phi):
    """Input phi in degrees, everything as numpy array"""
    phi = phi * np.pi / 180
    return cl * np.cos(phi) + cd * np.sin(phi)


t = {{0.22, 43.6}, {0.27, 38}, {0.3, 35}, {0.39, 29}, {0.45, 25}, {0.5, 22.5}, {0.6, 19.5}, {0.7, 17}, {0.8, 14.5}, {0.9, 13}, {0.95, 12.3}, {1, 12}}
