function y = my_fft(signal)
    commd = sprintf('./fft_calc %i exp.txt', length(signal));
    [~, ~] = system( commd );
    fileID = fopen('exp_fft.txt','r');
    y = fscanf(fileID,'%lf\n');

y = y(1:length(signal));