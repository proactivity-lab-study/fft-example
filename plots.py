import sys

import matplotlib.pyplot as plt
import numpy as np

# signal filename can be passed as: 
# python3 plots.py sin02
if len(sys.argv) > 1:
    file = sys.argv[1]
else:
    file = "sin01"

sampling_freq = 10000 # Hz
sampling_length = 6.5 #s

s = np.loadtxt(f'./signals/{file}.txt')
t = np.linspace(0, sampling_length, len(s))

fft_res = np.loadtxt(f'./signals/{file}_fft.txt', skiprows=1, usecols=[0,1])
nyq = round(len(fft_res)/2)

freq = fft_res[0:nyq, 0]
energy = fft_res[0:nyq, 1]

fig, (ax0, ax1) = plt.subplots(2, 1)
ax0.plot(t, s)
ax1.plot(freq, energy)

plt.show()