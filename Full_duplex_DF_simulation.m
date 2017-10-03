clear;close all;clc;j=1i;
%% Parameters
L = 30; % Number of Symbols
tau = 1; % Number of Delay
P_s = 25; % dB
P_r = 25; % dB
h_rr_gain = -20; % dB
h_sd_gain = -5; % dB
Eb_N0_all = 0:30; % dB
% Tx_number = 1e6; % Number of transmission
Tx_number = 17000; % Number of transmission
%% Main
BER_all = zeros(Tx_number,length(Eb_N0_all)); % Initialize
for s = 1:length(Eb_N0_all)
    Eb_N0 = Eb_N0_all(s);
    fprintf(2,['Eb_N0 = ', num2str(Eb_N0), '\n']); % Display
    h_sr = sqrt(0.5*(randn(Tx_number,1).^2+randn(Tx_number,1).^2)); % Rayleigh
    h_rr = sqrt(0.5*(randn(Tx_number,1).^2+randn(Tx_number,1).^2)); % Rayleigh
    Expected_value = sum(abs(h_sr).^2)/length(h_sr); % E[|h|^2]=1
    fprintf(1,['Expected_value = ', num2str(Expected_value), '\n']); % Display
    for Tx = 1:Tx_number
        %% Data_Payload generation
        M = 4; % QPSK
        Data_Payload = randi([0 M-1],L,1);
        %% Mapping
        x_s = pskmod(Data_Payload,M,pi/4);
        %% Rayleigh Fading Channel
        H_1 = sqrt(P_s)*h_sr(Tx)*[eye(L);zeros(tau,L)] + sqrt(P_r)*h_rr(Tx)*[zeros(tau,L);eye(L)];
        %% Noise
        sigma = sqrt(10^(-Eb_N0/10))/2;
        n_r = sigma*(randn(L+tau,1) + j*randn(L+tau,1));
        %% Received Signal
        y_r = H_1 * x_s + n_r;
        %% Zero-forcing Equalizer
        x_s_hat = H_1\y_r;
        %% DeMapping
        Data_Payload_hat = pskdemod(x_s_hat,M,pi/4);
        %% Error Calculation
        Error_number = sum(Data_Payload_hat ~= Data_Payload); % Number of Errors
        BER_all(Tx,s) = Error_number/L;
    end
end
BER = sum(BER_all)/Tx_number;
%% Graph
Eb_N0_Lin = 10.^(Eb_N0_all/10);
TheoryBER = 0.5.*(1-sqrt(Eb_N0_Lin./(Eb_N0_Lin+1)));
hold on;
semilogy(Eb_N0_all,TheoryBER,'bp-','LineWidth',2);
semilogy(Eb_N0_all,BER,'mx-','LineWidth',2);
hold off;
grid on;
legend('Rayleigh-Theory', 'Rayleigh-Simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER for QPSK modulation in Rayleigh channel');