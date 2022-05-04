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


def plot(mic, n=13):
    rpm = 8000

    # 51 channels
    # sampleRate: 51200
    # total measurements per sensor: 1024000
    # this means a experiment duration of 20 seconds
    spl2 = []
    PSD = np.zeros(5001)
    filename = ["measurements/noise/Bkg_U8_grid.h5", "measurements/noise/Motor_U0_rpm8000.h5",
                "measurements/noise/Propeller_U8_grid.h5"]
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

            data = measurements[n]
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

    # mic = [[69.16227277, 59.79007456, 38.8394522],
    #        [68.28505307, 59.78050416, 42.06984364],
    #        [72.77182762, 58.90286383, 45.26576453],
    #        [79.6853734, 62.25516431, 51.48410072],
    #        [84.26646802, 66.74432844, 57.47315295],
    #        [79.08249248, 64.21662571, 46.62326022],
    #        [73.57267351, 57.9725226, 43.89925439],
    #        [72.25592763, 56.12411415, 43.5055262],
    #        [73.77494749, 54.90854914, 41.72362697]]

    # micmean = np.mean(mic, axis=0)
    # plt.plot(1, mic[4][0], 'bo')
    # plt.plot(2, mic[4][1], 'bo')
    # plt.plot(3, mic[4][2], 'bo')
    #
    # plt.plot(1, micV[4][0], 'ro')
    # plt.plot(2, micV[4][1], 'ro')
    # plt.plot(3, micV[4][2], 'ro')

    plt.plot(1, mic[0], 'ro')
    plt.plot(2, mic[1], 'ro')
    plt.plot(3, mic[2], 'ro')

    plt.plot(fAdjusted, spl2)
    plt.xscale('log')
    # plt.ylim(ymin = -20, ymax = 80)
    # plt.xlim(xmin = 0.1, xmax = 100)
    plt.grid(which="both")
    plt.show()

    # Mic 25 Exact values- 65.3,47.23,55.7


micV = np.array([[68.27144740775242, 54.37963260949272, 38.274599443779216],
        [64.7214960814319, 54.22372493913146, 40.19270193961792],
        [68.83685572213328, 54.5371327346437, 41.05878253844993],
        [77.18757512326223, 62.08388021068666, 47.318978487235725],
        [82.08245238228133, 67.68624954216197, 53.23379805221197],
        [77.18757512326223, 62.08388021068666, 47.318978487235725],
        [68.83685572213328, 54.5371327346437, 41.05878253844993],
        [64.7214960814319, 54.22372493913146, 40.19270193961792],
        [68.27144740775242, 54.37963260949272, 38.274599443779216]])

micTN1 = np.array([[65.19796906398403, 52.59706506448667, 36.82324172270197],
          [64.4403627670812, 52.82919563161321, 38.899796318336854],
          [71.99968603732972, 56.39140435289853, 41.831916400286985],
          [78.61099286223039, 63.58807296475157, 48.77293167242275],
          [82.77952022074209, 68.45096398515241, 54.01530639139979],
          [78.61099286223039, 63.58807296475157, 48.77293167242275],
          [71.99968603732972, 56.39140435289853, 41.831916400286985],
          [64.4403627670812, 52.82919563161321, 38.899796318336854],
          [65.19796906398403, 52.59706506448667, 36.82324172270197]])

micTN2 = np.array([[70.6350207003376, 56.01251759022459, 39.58021104143517],
          [68.20102581326223, 56.310956596466454, 41.74580946394862],
          [65.89541721635096, 54.6307723140813, 41.68663378739783],
          [74.90504739587317, 60.027141345134666, 45.58044728070896],
          [80.63380001211193, 66.37380566479163, 52.03972260135459],
          [74.90504739587317, 60.027141345134666, 45.58044728070896],
          [65.89541721635096, 54.6307723140813, 41.68663378739783],
          [68.20102581326223, 56.310956596466454, 41.74580946394862],
          [70.6350207003376, 56.01251759022459, 39.58021104143517]])

micAve = np.log10((10**micTN2+10**micTN1+10**micV)/3)

plot(micAve[5])
