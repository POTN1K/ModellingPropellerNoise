import h5py


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


# Open files
filename = "Motor_U0_rpm8000.h5"
with h5py.File(filename, "r") as f:
    openKey(f, 0)
    print(readData(f, "/Acquisition/NumberOfChannels")[()])
    print(readData(f, "/Acquisition/NumberOfSamples")[()])
    print(readData(f, "/Acquisition/SampleRate")[()])
    print(readData(f, "/Acquisition/MicrophoneID")[()])
