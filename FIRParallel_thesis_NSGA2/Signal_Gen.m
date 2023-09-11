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

SNR_values = [1.1 2 5 10 20 30]  

for i = 1:numel(SNR_values)
  % Add a noise
  SNR = SNR_values(i);
  sine_noise = awgn(sine_wave,SNR, PWR = 'measured', seed=2);
  #sine_noise = sine_wave + Amp*sin(20*t) + Amp*sin(40*t) + Amp*sin(60*t);
  #sine_norm = sine_noise / max(abs(sine_noise));
  figure();plot(1:length(sine_noise), sine_noise);
  xlabel('\bf Time');
  ylabel('\bf Amplitude');
  title(sprintf("\b Sine wave + Noise: SNR of %d dB", SNR));
  
  ui = typecast(sine_noise, 'uint64');
  sine_noise_binary = dec2bin(ui,64)  
  yy = cellstr(sine_noise_binary);
  filename1 = sprintf("input_binary_%dSNR.data",SNR);
  fid = fopen(filename1, 'wt');
  fprintf(fid, '%8s \n', yy{:});
  disp('text file for signal finished');
  fclose(fid)
 
  filename2 = sprintf("input_%dSNR.data",SNR);
  csvwrite(filename2,sine_noise) 

end




