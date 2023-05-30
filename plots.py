import sys

import matplotlib.pyplot as plt
import numpy as np

# Signal filename can be passed as: 
# python3 plots.py sin02
if len(sys.argv) > 1:
    file = sys.argv[1]
else:
    file = "sin01"

sampling_freq = 10000 # Hz
sampling_length = 6.5 #s

s = np.loadtxt(f'./signals/{file}.txt')
t = np.linspace(0, sampling_length, len(s))

# Visualize original signal and fft of full signal 
fft_res = np.loadtxt(f'./signals/{file}_fft.txt', skiprows=1, usecols=[0,1])

x_lim = -1
x_axis_freqs = fft_res[0:x_lim, 0]
energy = fft_res[0:x_lim, 1]

fig, (ax0, ax1, ax2) = plt.subplots(3, 1)
ax0.set_title("Signal")
ax0.plot(t, s)
ax1.set_title("FFT on full signal", y=1.0, pad=-14)
ax1.plot(x_axis_freqs, energy)

# Visualize fft of windows
fft_windows_res = np.loadtxt(f'./signals/{file}_fft_only_windows.txt', skiprows=1)

# take mean of all windows for each frequency in x-axis
windows_mean = np.mean(fft_windows_res[0:-1, 1:-1], axis=1, dtype=np.float32)

x_axis_freqs = fft_windows_res[0:x_lim, 0]
energy = np.transpose(windows_mean)

ax2.set_title("FFT - moving windows (avg)", y=1.0, pad=-14)
ax2.plot(x_axis_freqs, energy)

plt.show()