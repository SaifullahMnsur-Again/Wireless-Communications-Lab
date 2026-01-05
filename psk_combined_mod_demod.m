clc;
clear all;
close;

M = [2, 4, 8, 16, 32];
numBits = 1e5;
SNR_db = 0:2:40;

Fs = 1e5;
maximumDopplerShift = 30;
pathDelays = 0;
averagePathGains = 0;
K_Factor = 10;

rayChan = comm.RayleighChannel(...
    'SampleRate', Fs, ...
    'PathDelays', pathDelays, ...
    'AveragePathGains', averagePathGains, ...
    'MaximumDopplerShift', maximumDopplerShift, ...
    'PathGainsOutputPort', true);

ricChan = comm.RicianChannel(...
    'SampleRate', Fs, ...
    'PathDelays', pathDelays, ...
    'AveragePathGains', averagePathGains, ...
    'MaximumDopplerShift', maximumDopplerShift, ...
    'KFactor', K_Factor, ...
    'PathGainsOutputPort', true);

for M_idx = M
    k = log2(M_idx);
    
    ber_awgn = [];
    ber_ray= [];
    ber_ric= [];
    
    data = randi([0, 1], numBits, 1);
    
    remainder = mod(numBits, k);
    if remainder ~= 0
        extra = k - remainder
        padding = zeros(extra, 1);
        data_dadded = [data; padding];
    else
        data_dadded = data;
    end
    
    symbol_count = length(data_dadded)/k;
    data_reshaped = reshape(data_dadded, k, symbol_count)';
    data_symbols = bi2de(data_reshaped, 'left-msb');
    
    for SNR = SNR_db
        modSignal = pskmod(data_symbols, M_idx);
        
        rx_awgn = awgn(modSignal, SNR, 'measured');
        demod_awgn = pskdemod(rx_awgn, M_idx);
        
        bits_awgn = de2bi(demod_awgn, k, 'left-msb')';
        bits_awgn = bits_awgn(:);
        [~, err_awgn] = biterr(data, bits_awgn(1:numBits));
        ber_awgn = [ber_awgn, err_awgn];
        
        
        reset(rayChan);
        [txFaded_ray, h_ray] = step(rayChan, modSignal);
        rx_ray = awgn(txFaded_ray, SNR, 'measured');
        
        rxEq_ray = rx_ray ./ squeeze(h_ray);
        
        demod_ray = pskdemod(rxEq_ray, M_idx);
        
        bits_ray = de2bi(demod_ray, k, 'left-msb')';
        bits_ray = bits_ray(:);
        [~, err_ray] = biterr(data, bits_ray(1:numBits));
        ber_ray = [ber_ray, err_ray];
        
        reset(ricChan);
        [txFaded_ric, h_ric] = step(ricChan, modSignal);
        rx_ric = awgn(txFaded_ric, SNR, 'measured');
        
        rxEq_ric = rx_ric./squeeze(h_ric);
        
        demod_ric = pskdemod(rxEq_ric, M_idx);
        
        bits_ric = de2bi(demod_ric, k, 'left-msb')';
        bits_ric = bits_ric(:);
        [~, err_ric] = biterr(data, bits_ric(1:numBits));
        ber_ric = [ber_ric, err_ric];
    end
    figure;
    semilogy(SNR_db, ber_awgn, 'o-', 'LineWidth', 2, 'DisplayName', 'AWGN');
    hold on;
    semilogy(SNR_db, ber_ray, 's-', 'LineWidth', 2, 'DisplayName', 'Rayleigh');
    semilogy(SNR_db, ber_ric, 'd-', 'LineWidth', 2, 'DisplayName', 'Rician');
    
    grid on;
    xlabel('SNR (db)');
    ylabel('Bit Error Rate (BER)');
    title(sprintf('Channel Comparison for %d-PSK', M_idx));
    legend('show');
end