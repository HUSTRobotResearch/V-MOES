import numpy as np
import os

# os.chdir("exp1")
map = "map1a"
'hypervolume' in dir()
allfile = os.listdir("C:/Users/ADMIN/Downloads/Multi/Multi/exp2/"+map)
# allfile.sort()
print(allfile)
box = []
res = [[], [], [], []]
i = 0

for file in allfile:
    f = open("C:/Users/ADMIN/Downloads/Multi/Multi/exp2/"+map+"/" + file, "r")

    flag = False

    while True:
        line = f.readline()
        tempstr = []
        while line.startswith("P") == True:
            flag = True
            mystr = line.split()
            for arr in mystr:
                if arr.replace('.', '', 1).isdigit():
                    tempstr.append(float(arr))
            res[i].append(tempstr)
            line = f.readline()
            tempstr = []
        if (line.startswith("-") == True):
            flag = False
        if ("" == line):
            print("file finished")
            break

    f.close()

    i += 1

ratio = [[], [], [], []]

for i in range(4):
    for j in range (4):
        count = 0
        for k in res[i]:
            for l in res[j]:
                if l[0] <= k[0] and l[1] <= k[1] and l[2] <= k[2]:
                    count += 1
                    break
        if len(res[i]) == 0:
            ratio[i].append(0.0)
        else:
            ratio[i].append(count / len(res[i]) * 100)

print("vmoes/\tmoes/\tmopso/\tnsgaii/")
print("vmoes", ratio[0])
print("moes", ratio[1])
print("mopso", ratio[2])
print("nsgaii", ratio[3])
