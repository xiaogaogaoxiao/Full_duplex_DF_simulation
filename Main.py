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
#Tx_number = 10**6 # Number of transmission
Tx_number = 17000 # Number of transmission
BER_all = np.zeros([np.size(Eb_N0_all),Tx_number]) # Initialize

# Main
def main():
    for s in range(0, np.size(Eb_N0_all), 1):
        Eb_N0 = Eb_N0_all[s]
        h_sr = Rayleigh_generate() # Rayleigh
        h_rr = Rayleigh_generate()*np.sqrt(h_rr_gain_Linear) # Rayleigh
        h_sd = Rayleigh_generate()*np.sqrt(h_sd_gain_Linear) # Rayleigh
        h_rd = Rayleigh_generate() # Rayleigh
        Expected_value_h_sr = Expected_value(h_sr) # E[|h|^2]=1
        Expected_value_h_rr = Expected_value(h_rr) # E[|h|^2]=0.01
        Expected_value_h_sd = Expected_value(h_sd) # E[|h|^2]=0.3162
        Expected_value_h_rd = Expected_value(h_rd) # E[|h|^2]=1

def Rayleigh_generate():
    return 1/np.sqrt(2)*((np.random.randn(1,Tx_number)+1j*np.random.randn(1,Tx_number)))

def Expected_value(x):
    return np.mean(np.abs(x)**2)

if __name__ == '__main__':
    cProfile.run('main()')

# End Program
end = time.time()
elapsed = end - start
print("Time taken: ", elapsed, "seconds.")