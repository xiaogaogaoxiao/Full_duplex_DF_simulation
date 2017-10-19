clear;close all;clc;j=1i;
tic;
parpool('local',8);
%% Parameters
L = 30; % Number of Symbols
tau = 1; % Number of Delay
P_s = 0; % dB
P_s_Linear = 10^(P_s/10); % Linear
P_r = 0; % dB
P_r_Linear = 10^(P_r/10); % Linear
h_rr_gain = -20; % dB
h_rr_gain_Linear = 10^(h_rr_gain/10);
h_sd_gain = -5; % dB
h_sd_gain_Linear = 10^(h_sd_gain/10);
Eb_N0_all = 0:2:30; % dB
Tx_number = 1e6; % Number of transmission
% Tx_number = 17000; % Number of transmission
%% Main
BER_all = zeros(Tx_number,length(Eb_N0_all)); % Initialize
for s = 1:length(Eb_N0_all) % SNR Loop
    Eb_N0 = Eb_N0_all(s);
    fprintf(2,['Eb_N0 = ', num2str(Eb_N0), '\n']); % Display
    h_sr = 1/sqrt(2)*((randn(Tx_number,1)+j*randn(Tx_number,1))); % Rayleigh
    h_rr = 1/sqrt(2)*((randn(Tx_number,1)+j*randn(Tx_number,1)))*sqrt(h_rr_gain_Linear); % Rayleigh
    h_sd = 1/sqrt(2)*((randn(Tx_number,1)+j*randn(Tx_number,1)))*sqrt(h_sd_gain_Linear); % Rayleigh
    h_rd = 1/sqrt(2)*((randn(Tx_number,1)+j*randn(Tx_number,1))); % Rayleigh
    Expected_value_h_sr = sum(abs(h_sr).^2)/length(h_sr); % E[|h|^2]=1
    Expected_value_h_rr = sum(abs(h_rr).^2)/length(h_rr); % E[|h|^2]=0.01
    Expected_value_h_sd = sum(abs(h_sd).^2)/length(h_sd); % E[|h|^2]=0.3162
    Expected_value_h_rd = sum(abs(h_rd).^2)/length(h_rd); % E[|h|^2]=1
    fprintf(1,['Expected_value = ', num2str(Expected_value_h_sr), '\n']); % Display
    parfor Tx = 1:Tx_number % Transmission Loop
        %% Data_Payload generation
        M = 4; % QPSK
        Data_Payload = randi([0 M-1],L,1);
        %% Mapping
        x_s = pskmod(Data_Payload,M,pi/4);
        %% Rayleigh Fading Channel (S-R)
        H_1 = sqrt(P_s_Linear)*h_sr(Tx)*[eye(L);zeros(tau,L)] + sqrt(P_r_Linear)*h_rr(Tx)*[zeros(tau,L);eye(L)];
        %% Noise (S-R)
        N0 = 10^(-Eb_N0/10);
        sigma = sqrt(N0/2);
        n_r = 1/sqrt(2)*sigma*(randn(L+tau,1) + j*randn(L+tau,1));
        %% Received Signal (S-R)
        y_r = H_1 * x_s + n_r;
        %% Zero-forcing Equalizer (S-R)
        H_1_hat = H_1;
        x_s_hat = H_1_hat\y_r; % inv(H_1_hat)*y_r
        %% Rayleigh Fading Channel (R-D)
        H_2 = sqrt(P_s_Linear)*h_sd(Tx)*[eye(L);zeros(tau,L)] + sqrt(P_r_Linear)*h_rd(Tx)*[zeros(tau,L);eye(L)];
        %% Noise (R-D)
        N0 = 10^(-Eb_N0/10);
        sigma = sqrt(N0/2);
        n_d = 1/sqrt(2)*sigma*(randn(L+tau,1) + j*randn(L+tau,1));
        %% Received Signal (R-D)
        y_d = sqrt(P_s_Linear)*h_sd(Tx)*[x_s;zeros(tau,1)] + sqrt(P_r_Linear)*h_rd(Tx)*[zeros(tau,1);x_s_hat] + n_d;
%         y_d = H_2 * x_s_hat + n_d;
        %% Zero-forcing Equalizer (R-D)
        H_2_hat = H_2;
        x_s_hat_hat = H_2_hat\y_d; % inv(H_2_hat)*y_d
        %% DeMapping
        Data_Payload_hat = pskdemod(x_s_hat_hat,M,pi/4);
        %% Error Calculation
        Error_number = sum(Data_Payload_hat ~= Data_Payload); % Number of Errors
        BER_all(Tx,s) = Error_number/L;
    end
end
BER = sum(BER_all)/Tx_number/2/2;
%% Graph
TheoryBERAWGN = 0.5*erfc(sqrt(10.^(Eb_N0_all/10))); % Theoretical AWGN BER
TheoryBER = 0.5.*(1-sqrt(10.^(Eb_N0_all/10)./(10.^(Eb_N0_all/10)+1))); % Theoretical Rayleigh BER
semilogy(Eb_N0_all,TheoryBERAWGN,'cd-','LineWidth',2);
hold on;
semilogy(Eb_N0_all,TheoryBER,'bp-','LineWidth',2);
semilogy(Eb_N0_all,BER,'mx-','LineWidth',2);
hold off;
grid on;
legend('AWGN','Rayleigh-Theory', 'Full-duplex-Simulation');
axis([0 30 10^-5 0.5]);
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER for QPSK modulation in Rayleigh channel');
%% Close Parpool
delete(gcp('nocreate'));
toc;