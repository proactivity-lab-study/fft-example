close all; clc ; clear;

fpath = './signals/smenete_1.txt';

fid = fopen(fpath, 'r');
signal = fscanf(fid, '%f');

%wd = hanning(length(signal));
%signal = signal .* wd;
    
disp("*** With Matlab fft implementation ***")
detector21(signal);

fprintf("\n\n*** With custom fft implementation ***\n")
detector_custom_fft(signal);

%%
close all;
fpath = './signals/smenete_2.txt';

fid = fopen(fpath, 'r');
signal = fscanf(fid, '%f');


%wd = parzenwin(length(signal));
%signal_wd = signal .* wd;

disp("*** With Matlab fft implementation ***")
detector21(signal, 3000, 0.55, 1);

%detector21(signal, 3000, 0.65, 1.5);

fprintf("\n\n*** With custom fft implementation ***\n")
detector_custom_fft(signal, 3000, 0.55, 0.5);

%detector_custom_fft(signal, 3000, 0.6, 1);


%%
close all; clc ; clear;

% sig_191108_163923.txt';
 
%sdir = '/home/linomp/code/fft-example/Andmed/08112019/microwave/MICROWAVE';
%sdir = '/home/linomp/code/fft-example/Andmed/09112019/microwave/MICROWAVE';
sdir = '/home/linomp/code/fft-example/Andmed/30102019/microwave/MICROWAVE';

files = dir(sprintf('%s/*.txt', sdir));

for i=1:length(files)
    try
        file = erase(files(i).name, '.txt');

        fpath = sprintf('%s/%s.txt', sdir, file);
        fid = fopen(fpath, 'r');
        signal = fscanf(fid, '%f');

        disp("*** With Matlab fft implementation ***")
        [~,~, ~, fig] = detector21(signal, 3000, 0.5, 2);
        saveas(fig, sprintf('30102019/%s.png',file));

        fprintf("\n\n*** With custom fft implementation ***\n")
        [~,~, ~, fig] = detector_custom_fft(signal, 3000, 0.5, 0.5);
        saveas(fig, sprintf('30102019/%s_custom_fft.png',file));
         close all;
    catch
        close all;
    end
end

