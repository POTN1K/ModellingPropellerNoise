"""In this file we compare the experimental and Hanson Results"""

import h5py
import matplotlib.pyplot as plt
import scipy
from scipy import signal
import math
import numpy as np


def openKey(data, key):
    if (key == 0):
        newKeys = list(data.keys())
        # print("key: " + str(key) + " childKeys:" + str(newKeys))
        print("reading main file, contents: " + str(newKeys))
        for x in range(0, len(newKeys)):
            # print(newKeys[x])
            openKey(data, newKeys[x])
    else:

        newData = data[key]

        if hasattr(newData, 'keys'):
            newKeys = list(newData.keys())
            # print("key: " + str(key) + " childKeys:" + str(newKeys))
            print("reading folder '" + str(key) + "', contents: " + str(newKeys))
            for x in range(0, len(newKeys)):
                # print(newKeys[x])
                # newparent = str(parent) + "/" +  str(key)
                openKey(newData, newKeys[x])
        else:
            print(str(newData.name))
            # print(newData)
    return 0


def readData(data, path):
    pathArr = path.split("/")
    tempData = data
    if (pathArr[0] == ""):
        pathArr.pop(0)
    # print(pathArr)
    for x in range(0, len(pathArr)):
        newData = tempData[pathArr[x]]
        # print(newData)
        tempData = newData
    # print(tempData[()])
    return tempData


def numInt(x, y):
    sum = 0
    length = 0
    for i in range(0, len(x) - 1):
        df = x[i + 1] - x[i]
        v = y[i] * df
        sum += v
        length += x[i]
    sum = sum / length
    return sum


rpm = 8000

# 51 channels
# sampleRate: 51200
# total measurements per sensor: 1024000
# this means a experiment duration of 20 seconds
spl2 = []
PSD = np.zeros(5001)
filename = ["measurements/noise/Bkg_U8_nogrid.h5", "measurements/noise/Motor_U0_rpm8000.h5",
            "measurements/noise/Propeller_U8_nogrid.h5"]
# filename = "Bkg_U8_grid.h5"
# filename = "Bkg_U8_nogrid.h5"
# filename="Motor_U0_rpm8000.h5"
# filename = 'Bkg_U8_grid.h5'
for filename in filename:
    with h5py.File(filename, "r") as f:
        openKey(f, 0)
        # print(readData(f, "/Acquisition/NumberOfChannels")[()])
        # print(readData(f, "/Acquisition/NumberOfSamples")[()])
        # print(readData(f, "/Acquisition/SampleRate")[()])
        measurements = readData(f, "/AcousticData/Measurement")[()]

        fs = 51200
        window = 10000
        overlap = window / 2
        nfft = window
        p_ref = 20 * 10 ** -6
        df = fs / window
        '''f_start = 20'''
        f_start = 20
        f_end = 200

        data = measurements[25]
        PSD1 = PSD
        # print(PSD1)
        f, PSD = scipy.signal.welch(data, fs=fs, window='hamming', nperseg=window, nfft=nfft, noverlap=overlap)
        bpf = (rpm / 60) * 2

        PSD2 = (PSD - PSD1)
        spl = []
        #     print(PSD)
        #    print(PSD2)
        for pressureLevel in PSD:
            valueDecibels = 10 * math.log10(pressureLevel / (p_ref ** 2))  # +10 shifte wrt Luc van Beel
            spl.append(valueDecibels)
        #        spl = np.array(spl)
        #       spl2 = np.array(spl2)
        fAdjusted = []
        for freq in f:
            fAdjusted.append(freq / bpf)

        fAdjusted = np.array(fAdjusted)
        #  fAdjusted [fAdjusted > 0.1]
        OASPL = 10 * math.log10(np.std(data) / (p_ref))
    # print(f)
    # print(OASPL)
    # print(numInt(f,spl))
    # spl2 = np.array(spl) - np.array(spl2)

    plt.plot(fAdjusted, spl)

    # print(np.size(PSD2))
# print(PSD2)
# print(spl2)
# plt.plot(fAdjusted, spl2)
for i in np.where(np.array(PSD2) < 10 ** -8):
    PSD2[i] = 20 * 10 ** -6
for pressureLevel in PSD2:
    valueDecibels = 10 * math.log10(pressureLevel / (p_ref ** 2))  # +10 shifte wrt Luc van Beel
    spl2.append(valueDecibels)

mic = [[95.32667909, 74.68890744, 60.70080301],
       [97.49193266, 78.33964198, 65.41170831],
       [99.30079103, 81.38312807, 69.37623245],
       [100.52393853, 83.43759727, 72.06793225],
       [100.95922676, 84.16800581, 73.02749453],
       [99.07462501, 86.98971912, 68.08994655],
       [97.83321272, 85.0216877, 65.23549884],
       [95.9965093, 82.11393138, 61.00411331],
       [93.79653672, 78.63941566, 55.92497085]]

# mic = np.array([[66.9374063885016, 52.71483119068122, 36.55945810157958],
#                 [69.14969071197915, 56.37397810382691, 41.65302713946032],
#                 [70.99692474217275, 59.42529927423304, 45.89623124996426],
#                 [72.24563635965032, 61.48556948283519, 48.758727071510776],
#                 [72.68995884301663, 62.21815290684128, 49.77601820976848],
#                 [72.24563635965032, 61.48556948283519, 48.758727071510776],
#                 [70.99692474217275, 59.42529927423304, 45.89623124996426],
#                 [69.14969071197915, 56.37397810382691, 41.65302713946032],
#                 [66.9374063885016, 52.71483119068122, 36.55945810157958]])

micmean = np.mean(mic, axis=0)
plt.plot(1, mic[0][0], 'bo')
plt.plot(2, mic[0][1], 'bo')
plt.plot(3, mic[0][2], 'bo')

plt.plot(fAdjusted, spl2)
plt.xscale('log')
# plt.ylim(ymin = -20, ymax = 80)
# plt.xlim(xmin = 0.1, xmax = 100)
plt.grid(which="both")
plt.show()
