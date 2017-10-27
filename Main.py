import numpy as np
import matplotlib.pyplot as plt
import cProfile
import time

start = time.time()

# Parameters
L = 30 # Number of Symbols
tau = 1 # Number of Delay
P_s = 0 # dB
P_s_Linear = 10**(P_s/10) # Linear
P_r = 0 # dB
P_r_Linear = 10**(P_r/10) # Linear
h_rr_gain = -20 # dB
h_rr_gain_Linear = 10**(h_rr_gain/10)
h_sd_gain = -5 # dB
h_sd_gain_Linear = 10**(h_sd_gain/10)
Eb_N0_all =  np.arange(0,32,2) # dB
Tx_number = 10**6 # Number of transmission
BER_all = np.zeros([np.size(Eb_N0_all),Tx_number]) # Initialize

# Main
for s in range(0, np.size(Eb_N0_all), 1): # SNR Loop
    Eb_N0 = Eb_N0_all[s]
    h_sr = 1/np.sqrt(2)*((np.random.randn(1,Tx_number)+1j*np.random.randn(1,Tx_number))) # Rayleigh
    h_rr = 1/np.sqrt(2)*((np.random.randn(1,Tx_number)+1j*np.random.randn(1,Tx_number)))*np.sqrt(h_rr_gain_Linear) # Rayleigh
    h_sd = 1/np.sqrt(2)*((np.random.randn(1,Tx_number)+1j*np.random.randn(1,Tx_number)))*np.sqrt(h_sd_gain_Linear) # Rayleigh
    h_rd = 1/np.sqrt(2)*((np.random.randn(1,Tx_number)+1j*np.random.randn(1,Tx_number))) # Rayleigh
    Expected_value_h_sr = np.mean(np.abs(h_sr)**2) # E[|h|^2]=1
    Expected_value_h_rr = np.mean(np.abs(h_rr)**2) # E[|h|^2]=0.01
    Expected_value_h_sd = np.mean(np.abs(h_sd)**2) # E[|h|^2]=0.3162
    Expected_value_h_rd = np.mean(np.abs(h_rd)**2) # E[|h|^2]=1

    for Tx in range(0, Tx_number, 1): # Transmission Loop
        # Data_Payload generation
        M = 4 # QPSK
        Data_Payload = np.random.randint(M, size=L)

# End Program
end = time.time()
elapsed = end - start
print("Time taken: ", elapsed, "seconds.")