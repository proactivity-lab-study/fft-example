close all; clc ; clear;

fid = fopen('./signals/smenete_1.txt', 'r');
signal = fscanf(fid, '%f');


disp("*** With Matlab fft implementation ***")
detector21(signal);

fprintf("\n\n*** With custom fft implementation ***\n")
detector_custom_fft(signal);

%%
fid = fopen('./signals/smenete_2.txt', 'r');
signal = fscanf(fid, '%f');


disp("*** With Matlab fft implementation ***")
detector21(signal);

fprintf("\n\n*** With custom fft implementation ***\n")
detector_custom_fft(signal, 3000, 0.35, 1);
