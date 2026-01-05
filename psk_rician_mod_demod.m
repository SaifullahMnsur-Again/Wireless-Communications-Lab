clc;
clear all;
close;

M = [2, 4, 8, 16, 32];
numBits = 1e5;

SNR_db = 0:2:40;

maxDopplerShift = 30;
Fs = 1e5;
pathDelay = 0;
averagePathGain = 0;
K_Factor = 10;

ricChan = comm.RicianChannel(...
    'SampleRate', Fs, ...
    'PathDelays', pathDelay, ...
    'AveragePathGains', averagePathGain, ...
    'MaximumDopplerShift', maxDopplerShift, ...
    'KFactor', K_Factor, ...
    'PathGainsOutputPort', true);

figure();

for M_idx = M
    ber_vals = [];
    k = log2(M_idx);
    
    data = randi([0, 1], numBits, 1);
    
    remainder = mod(numBits, k);
    
    if remainder ~= 0
        extra = k - remainder;
        padding = zeros(extra, 1);
        dataPadded = [data; padding];
    else
        dataPadded = data;
    end
    
    symbolCount = length(dataPadded)/k;
    reshapedData = reshape(dataPadded, k, symbolCount)';
    dataSymbols = bi2de(reshapedData, 'left-msb');
    
    for SNR = SNR_db
        % Correct Order of Blocks: Modulation -> Rayleigh Fading -> AWGN -> Demodulation.
        
        % block-1 modulation
        modSignal = pskmod(dataSymbols, M_idx);
        
        reset(ricChan); % reset for new snr loop
        
        % block-2 Rayleigh Fading
        [txFaded, pathGains] = step(ricChan, modSignal); % step applies rayleigh channel on modded signal
        
        % block-3 AWGN
        rxNoisy = awgn(txFaded, SNR, 'measured');
        
        h = squeeze(pathGains); % remove extra dimensions
        rxEqualized = rxNoisy ./ h; % we must divide by the channel path gains (h) to fix phase rotation
        
        % block-4 Demodulation
        demodSignal = pskdemod(rxEqualized, M_idx);
        
        demodBits = de2bi(demodSignal, k, 'left-msb')';
        demodBitsFlat = demodBits(:);
        demodBitsWithoutPadding = demodBitsFlat(1:numBits);
        
        [~, ber] = biterr(data, demodBitsWithoutPadding);
        ber_vals = [ber_vals, ber];
    end
    
    semilogy(SNR_db, ber_vals, 'o-', 'LineWidth', 2);
    hold on;
end

grid on;
xlabel('SNR (db)');
ylabel('Bit Error Rate (BER)');
legend('BPSK', 'QPSK', '8-PSK', '16-PSK', '32-PSK')
title('BER vs SNR for Different PSK Modulation (Rician)');