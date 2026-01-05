clc;
clear all
close;

M = [2, 4, 8, 16];
numBits = 1e4;
SNR_db = 0:2:30;

Fs = 1e5;
nsamp = 16; % samples per symbol

% Symbol Rate = Fs/nsamp = 6250 Hz.
% Frequency separation must be > Symbol Rate for clean orthogonality.
freq_sep = Fs/nsamp + 1; % frequency seperator

maximumDopplerShift = 30;
pathDelays = 0;
averagePathGains = 0;
K_Factor = 10;

ricChan = comm.RicianChannel(...
    'SampleRate', Fs,...
    'PathDelays', pathDelays,...
    'AveragePathGains', averagePathGains,...
    'MaximumDopplerShift', maximumDopplerShift,...
    'KFactor', K_Factor, ...
    'PathGainsOutputPort', true);

figure()

for M_idx = M
    ber_vals = []
    k = log2(M_idx);
    
    data = randi([0, 1], numBits, 1);
    
    remainder = mod(numBits, k);
    if remainder ~= 0
        extra = k - remainder;
        padding = zeros(extra, 1);
        data_padded = [data; padding];
    else
        data_padded = data;
    end
    symbol_count = length(data_padded)/k;
    data_reshaped = reshape(data_padded, k, symbol_count)';
    data_symbols = bi2de(data_reshaped, 'left-msb');
    
    for SNR = SNR_db
        mod_signal = fskmod(data_symbols, M_idx, freq_sep, nsamp, Fs);
        
        reset(ricChan);
        [txFaded_signal, h_ric] = step(ricChan, mod_signal);
        rx_ric = awgn(txFaded_signal, SNR, 'measured');
        
        rxEq_symbols = rx_ric ./ squeeze(h_ric);
        
        demod_symbols = fskdemod(rx_symbols, M_idx, freq_sep, nsamp, Fs);
        
        demod_bits = de2bi(demod_symbols, k, 'left-msb')';
        demod_bits = demod_bits(:);
        [~, ber] = biterr(data, demod_bits(1:numBits));
        ber_vals = [ber_vals, ber];
    end
    semilogy(SNR_db, ber_vals, 'o-', 'LineWidth', 2);
    hold on;
end

grid on;
xlabel('SNR (db)');
ylabel('BER');
legend('2-FSK', '4-FSK', '8-FSK', '16-FSK');
title('BER vs SNR for FSK (Rician)');
