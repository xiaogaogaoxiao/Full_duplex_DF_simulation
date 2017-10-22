import numpy as np
import matplotlib.pyplot as plt
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
Eb_N0_all =  np.arange(0,31,2) # dB
Tx_number = 10**6 # Number of transmission
# Tx_number = 17000 # Number of transmission

# Main
BER_all = np.zeros([Tx_number,np.size(Eb_N0_all)]) # Initialize

# End Program
end = time.time()
elapsed = end - start
print("Time taken: ", elapsed, "seconds.")