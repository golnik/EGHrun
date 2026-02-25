import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys
import operator
import re


## norm
norm=1.1504905808866939
print(norm)


## PES scan
s = """
Calculation is performed. Results: [[-265.2799031477, -264.8812155052, -264.8393039783], 1, [0, 1]]
Calculation is performed. Results: [[-265.2819650567, -264.8795094292, -264.8421863181], 3, [0, 0]]
Calculation is performed. Results: [[-265.2713753772, -264.8800777188, -264.8309850111], 3, [0, 2]]
Calculation is performed. Results: [[-265.2799031477, -264.8812155052, -264.8393039783], 3, [0, 1]]
Calculation is performed. Results: [[-265.2819650567, -264.8795094292, -264.8421863181], 0]
Calculation is performed. Results: [[-265.2722266911, -264.8743208221, -264.837153478], 3, [0, -2]]
Calculation is performed. Results: [[-265.2506751089, -264.8696437005, -264.8111505408], 3, [0, 3]]
Calculation is performed. Results: [[-265.252139879, -264.8623809302, -264.8196427317], 3, [0, -3]]
Calculation is performed. Results: [[-265.2802569117, -264.8780397981, -264.8424748792], 1, [0, -1]]
Calculation is performed. Results: [[-265.2802569117, -264.8780397981, -264.8424748792], 3, [0, -1]]


"""
    
# ~ print(s)
x = re.findall("([0-9]+\.[0-9]+|-[0-9]+\.[0-9]+)",s)
print(x)


GS=[]
I1=[]
I2=[]
for i in range(int(len(x)/3)):
    GS.append(float(x[i*3]))
    I1.append(float(x[i*3+1]))
    I2.append(float(x[i*3+2]))
    
x = re.findall("(\[[0-9], [0-9]\]|\[[0-9], -[0-9]\])",s)
px=[]
for i in range(len(x)):
    if i==4:
        px.append(0)
    f = re.findall("[0-9]|-[0-9]",x[i])
    px.append( float(f[1])*norm )

px = np.array(px)

## create subplots
fig, ax = plt.subplots(2,1)
fig.tight_layout(pad=3.0)
plt.subplots_adjust(left=0.2)

## do polynomial fits
x = np.linspace( np.min(px), np.max(px) )
f_GS = np.polyfit(px, GS, 2)
f_I1 = np.polyfit(px, I1, 2)
f_I2 = np.polyfit(px, I2, 2)


## GS
ax = plt.subplot(2,1,2);
plt.scatter(px,GS)
plt.plot(x,f_GS[0]*x**2+f_GS[1]*x+f_GS[2])
ax.title.set_text("Ground state")
ax.set_xlabel("dr (A)")
ax.set_ylabel("E (a.u.)")
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.ylim([np.min(GS)-.01,np.min(GS)+.07])


## I1
ax = plt.subplot(2,1,1);
ax.title.set_text("Ionic states")
plt.plot(x,f_I1[0]*x**2+f_I1[1]*x+f_I1[2])
plt.scatter(px,I1)
ax.set_xlabel("dr (A)")
ax.set_ylabel("E (a.u.)")
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.3f'))


## I2 
plt.plot(x,f_I2[0]*x**2+f_I2[1]*x+f_I2[2])
plt.scatter(px,I2)
ax.set_xlabel("dr (A)")
ax.set_ylabel("E (a.u.)")
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

plt.savefig("allPES.png")
