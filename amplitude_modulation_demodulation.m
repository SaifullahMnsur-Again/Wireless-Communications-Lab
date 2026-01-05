clc;
close all;
clear;


Fs = 100e3; % sampling frequency
dt = 1/Fs; % time difference
t = 0 : dt : 0.005; % time

% message signal
Am = 5;
Fm = 1000;
mt = Am * sin(2*pi*Fm*t);

% carrier signal
Ac = 10;
Fc = 10000;
ct = Ac * sin(2*pi*Fc*t);

% modulated signal
st = (Ac + mt) .* sin(2*pi*Fc*t);

% demodulated signal
demodulated_st = (Am * Ac .* sin(2*pi*Fm*t)) / Ac;

figure(1)
subplot(4, 1, 1);
plot(t, mt);
title('Message signal');
subplot(4, 1, 2);
plot(t, ct);
title('Carrier signal');
subplot(4, 1, 3);
plot(t, st);
title('AM Modulated signal');
subplot(4, 1, 4);
plot(t, demodulated_st);
title('AM Demodulated message signal');

