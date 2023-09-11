% Generate a sine wave
close all; clear all;
pkg load communications
randn("state", "reset")
randn("seed", 2)

fs =10.15;
Amp = 1;
t = 0:1/fs:4*pi; % time vector

sine_wave = Amp*sin(t);
figure();
plot(t, sine_wave);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Sine wave');




% Add a noise
SNR = 5;
#sine_noise = awgn(sine_wave,SNR, PWR = 'measured', seed=2)
sine_noise = sine_wave + Amp*sin(20*t) + Amp*sin(40*t) + Amp*sin(60*t);
sine_norm = sine_noise / max(abs(sine_noise));
figure();plot(1:length(sine_norm), sine_norm);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title(sprintf("\b Sine wave + Noise: SNR of %d dB", SNR));

% Convert from real to integers
total_wordlength = 64;
scaling = 7;
sine_noise_integers = round(sine_norm.*(2^scaling));
figure();plot(1:length(sine_noise_integers), sine_noise_integers);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Sine wave + Noise : Scaled Signal');

% Convert from integers to binary
sine_noise_in_binary_signed = dec2bin(mod((sine_noise_integers),2^total_wordlength),total_wordlength);
#sine_noise_in_binary_signed = dec2bin(sine_noise_integers,64);
yy = cellstr(sine_noise_in_binary_signed);
fid = fopen('signal_noisy_128_sine.data', 'wt');
fprintf(fid, '%8s \n', yy{:});
disp('text file for signal finished');
fclose(fid)

% Scaled original sine wave
target = round(sine_wave.*(2^scaling));
figure();plot(1:length(target), target);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Sine wave: Scaled Target Signal');

csvwrite('signal_noisy_128_sine_noscale.data',sine_noise) 
csvwrite('target_128_sine_noscale.data',sine_wave)
csvwrite('target_scale.data',target)
csvwrite('input_data.data',sine_noise_integers) 


foo = typecast( uint16(bin2dec(sine_noise_in_binary_signed)), 'int16')

