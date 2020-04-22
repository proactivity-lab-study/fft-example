# FFT-Example

Example of usage of fast Fourier Transform algorithm.

FFT algorithm from http://paulbourke.net/miscellaneous/dft/ (2020 April 22)
The algorithm expects input signal length value to be a power of two.
The algorithm performs calculations in place (ie on the buffer given to it).

Test signals are in directory signals. There are four signals with length 6.5s
(sin01, sin02, sin03, random) and two signals with length 4s (sin800, sin2367).

Sampling rate for all test signals is 10 kHz. All signals are normalized -1...1
except sin03 which is -1.1...1.1. 
 * sin01 - pure sine wave 800 Hz 4s long from beginning, then 2.5s of silence
 * sin02 - pure sine wave 2367 Hz 4s long at the end, beginning 2.5s is silence
 * sin03 - sum of signals sin01, sin02, random
 * random - 6.5s of white noise (or close to wn)
 * sin800 - 4s sine wave 800 Hz
 * sin2367 - 4s sine wave 2367 Hz

# Build
g++ -o fft_calc main.c

# Usage
./fft_calc signal/sin01 signal/sin02
./fft_calc signal/sin03
