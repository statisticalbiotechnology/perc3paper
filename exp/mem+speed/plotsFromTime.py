#!/usr/env python3

import glob
import re
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

sns.set_context("poster")

resDir = "./res/"

runtimes = {}
for filename in glob.iglob(f'{resDir }/*.time'):
    match = re.search('.+res.per.([0-9]+?).([0-9]+?)k.+', filename)
    if match:
        ver = match.group(1)
        ver = ver[:1] + "." + ver[1:]
        size = int(match.group(2))*1000
        if not ver in runtimes.keys():
            runtimes[ver] = {}
        if not size in runtimes[ver].keys():
            runtimes[ver][size] = {
                'user':float("inf"),
                'system':float("inf"),
                'wall':float("inf"),
                'cpu':float("inf"),
                'resident':float("inf"),            
            }
        # print(ver,size)
        with open(filename,"r") as file:
            first=file.readline()
            match = re.search('([0-9.+-EE]+)user ([0-9.+-EE]+)system ([0-9]+:)?([0-9]+):([0-9.+-eE]+)elapsed (([0-9.]+))%CPU \([0-9]*avgtext.[0-9]*avgdata ([0-9]+)maxresident\)k', first)
            if match:
                user = float(match.group(1))
                system = float(match.group(2))
                #print(filename,first)
                wall = int(match.group(4))*60. + float(match.group(5))
                if match.group(3):
                    wall += int(match.group(3)[:-1])*3600.
                cpu = int(match.group(6))
                resident = float(match.group(8))/(1024.*1024.)
                #print(filename,resident)
                runtimes[ver][size]['user'] = min(runtimes[ver][size]['user'], user)
                runtimes[ver][size]['system'] = min(runtimes[ver][size]['system'], system)
                runtimes[ver][size]['wall'] = min(runtimes[ver][size]['wall'], wall)
                # runtimes[ver][size]['elapsed'] = min(runtimes[ver][size]['elapsed'], time.mktime(elapsed))
                runtimes[ver][size]['cpu'] = min(runtimes[ver][size]['cpu'], cpu)
                runtimes[ver][size]['resident'] = min(runtimes[ver][size]['resident'], resident)
# print (runtimes)
vers = [ver for ver in runtimes.keys() for size in runtimes["2.8"].keys() ]
sizes = [size for ver in runtimes.keys() for size in runtimes["2.8"].keys() ]
time = pd.DataFrame({
    "Version": vers,
    "size": sizes,
    "user": [runtimes[ver][size]['user'] for ver in runtimes.keys() for size in runtimes["2.8"].keys()],
    "system": [runtimes[ver][size]['system'] for ver in runtimes.keys() for size in runtimes["2.8"].keys()],
    "wall": [runtimes[ver][size]['wall'] for ver in runtimes.keys() for size in runtimes["2.8"].keys()],
    "cpu": [runtimes[ver][size]['cpu'] for ver in runtimes.keys() for size in runtimes["2.8"].keys()],
    "resident": [runtimes[ver][size]['resident'] for ver in runtimes.keys() for size in runtimes["2.8"].keys()]
})
time["cpu_rate"]=time["size"]/(time["user"]+time["system"])
time["wall_rate"]=time["size"]/time["wall"]

f, ax = plt.subplots(figsize=(10, 10))
#ax.set(xscale="log", yscale="log")
sns.lineplot(data=time, x="size", y="user", hue = "Version", 
    hue_order = ['2.8', '3.0', '3.3', '3.6'], marker="o", ax=ax)
plt.xlabel("Number of PSMs")
plt.ylabel("User time [s]")
plt.savefig("img/user.png")

f, ax = plt.subplots(figsize=(10, 10))
sns.lineplot(data=time, x="size", y="system", hue = "Version", 
    hue_order = ['3.6', '3.3', '3.0', '2.8'], marker="o", ax=ax)
plt.xlabel("Number of PSMs")
plt.ylabel("System time [s]")
plt.savefig("img/system.png")
f, ax = plt.subplots(figsize=(10, 10))

sns.lineplot(data=time, x="size", y="wall", hue = "Version", 
    hue_order = ['2.8', '3.0', '3.3', '3.6'], marker="o", ax=ax)
plt.xlabel("Number of PSMs")
plt.ylabel("Wall time [s]")
plt.savefig("img/wall.png")

f, ax = plt.subplots(figsize=(10, 10))
sns.lineplot(data=time, x="size", y="resident", hue = "Version", 
    hue_order = ['2.8', '3.0', '3.3', '3.6'], marker="o", ax=ax)
plt.xlabel("Number of PSMs")
plt.ylabel("Peak Resident Memory Size [GB]")
plt.savefig("img/memory.png")

f, ax = plt.subplots(figsize=(10, 10))
sns.lineplot(data=time[time["size"]>5000], x="size", y="cpu_rate", hue = "Version", 
    hue_order = ['3.6', '3.3', '3.0', '2.8'], marker="o", ax=ax)
plt.xlabel("Number of PSMs")    
plt.ylabel("PSMs Processed per CPU second")
plt.ylim(0,5000)
plt.savefig("img/rate.png")

f, ax = plt.subplots(figsize=(10, 10))
sns.lineplot(data=time[time["size"]>5000], x="size", y="wall_rate", hue = "Version", 
    hue_order = ['3.6', '3.3', '3.0', '2.8'], marker="o", ax=ax)
plt.xlabel("Number of PSMs")
plt.ylabel("PSMs Processed per Wall second")
#plt.ylim(0,5000)
plt.savefig("img/wallrate.png")
