clc;
clearvars;
close all;

[mic_1, Fs_1] = audioread("Corrupted_Speech.wav");
[mic_2, Fs_2] = audioread("White_Noise.wav");
[clean, Fs_3] = audioread("cleanspeech.wav");

win_size = 2^8;
adapt = 0.0026;
nFrames = fix(min(length(mic_1), length(mic_2))/win_size);
initialFreqResp = zeros(win_size, 1);

% Fixing length of reduced noise signal same as 'cleanspeech' signal
opSig_f = zeros(length(clean), 1);
opSig_t = zeros(length(clean), 1);

for k = 1 : nFrames
    % Compute indices for current frame
    j = (1:win_size) + (win_size*(k-1));
    D = fft(mic_1(j), win_size);
    X = fft(mic_2(j), win_size);
    X_diag = diag(X); 
   
    if k == 1
        % filtered output signal
        oS_f(1:win_size, k) = D - (X_diag*initialFreqResp);
        oS_t(1:win_size, k) = ifft(oS_f(1:win_size, k));
        
        % We now again calculate reduced noise speech signal as a column-vector as that would help with the SNR calculation
    
        opSig_f(j, 1) = D - (X_diag*initialFreqResp);
        opSig_t(j, 1) = ifft(opSig_f(j, 1));

        % frequency response of the filter
        B(:, k) = initialFreqResp + (2*adapt*X_diag'*oS_f(1:win_size, k));

    else
        oS_f(1:win_size, k) = D - (X_diag*B(:, k-1));
        oS_t(1:win_size, k) = ifft(oS_f(1:win_size, k));
    
        opSig_f(j, 1) = D - (X_diag*B(:, k-1));
        opSig_t(j, 1) = ifft(opSig_f(j, 1));

        B(:, k) = B(:, k-1) + (2*adapt*X_diag'*oS_f(1:win_size, k));
    end

end

SNR_before = 10 * log10((clean'*clean)/((clean-mic_1)'*(clean-mic_1)));
SNR_after = 10 * log10((clean'*clean)/((clean-opSig_t)'*(clean-opSig_t)));
SNR_change = SNR_after - SNR_before;

for p = 1:nFrames
    % Calculating energy of each frame to find energy convergence
    E(1, p) = 10 * log10(oS_f(:, p)'*oS_f(:, p));

end

figure;
plot(E);
msg = sprintf("CONVERGENCE CURVE (N=%d, U=%f)", win_size, adapt);
title(msg);
xlabel("Frame number");
ylabel("Energy (dB)");

audiowrite("Filtered_Speech.wav", opSig_t, Fs_1);