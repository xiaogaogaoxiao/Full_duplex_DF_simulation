clear;close all;clc;j=1i;
%% Parameters
M = 4; % QPSK
L = 30;
tau = 1;
P_s = 25; % dB
P_r = 25; % dB
h_rr = -20; % dB
h_sd = -5; % dB
%% Main
data_Payload = randi([0 M-1],1,L);
s = pskmod(data_Payload,M,pi/4);
data_Payload_hat = pskdemod(s,M,pi/4);
plot(s,'x');