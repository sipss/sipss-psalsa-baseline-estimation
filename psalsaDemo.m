clear
close all
clc

%load a data set with a synthetic IMS spectra which includes an artificial baseline
load SynSpec2.mat % Load the dataset: variables inside are SynSpec_bl and bl
                   %bl is a synthetic sinusoid baseline
                   %SynSpec_bl is a synthetic noisy IMS spectra with the bl systhetic baseline added
figure
plot(SynSpec_bl), xlabel('samples'), ylabel('Spectrum Amplitude (a.u.)') %plot the IMS spectrum with baseline

%Apply psalsa for baseline correction using default values
[Xc, Z1] = psalsa(SynSpec_bl); %Generates the baseline Z and the corrected spectrum Xc
%plot the corrected spectrum
figure
plot(Xc), xlabel('samples'), ylabel('Spectrum Amplitude (a.u.)')

%Apply psalsa for baseline correction with a different Lambda parameter
[Xc,Z2] = psalsa(SynSpec_bl,'Lambda',1e3); %Changing the smoothing parameter
%plot the corrected spectrum
figure
plot(Xc), xlabel('samples'), ylabel('Spectrum Amplitude (a.u.)')

figure
title('baseline estimation comparative')
plot(bl) %plot the synthetic baseline
hold on
plot(Z1) %plot the estimated baseline with the default parameters
plot(Z2) %plot the estimated baseline with the smoothing Lambda parameter set to 1e3
xlabel('samples'), ylabel('Amplitude a.u.')
legend ('Synthetic baselie','Lambda = 1e6','Lambda = 1e3')

