from core import *
import time
import pickle

C_ = 121e-12
l_phys = 0.5
freq = 10e9
q = 0.190

filename = "last_simulation.pkl"

# Get target values from MATLAB
Idc_A = np.array([-0.03325, -0.030875, -0.0285, -0.026125, -0.02375, -0.021375, -0.019, -0.016625, -0.01425, -0.011875, -0.0095, -0.007125, -0.00475, -0.002375, 0, 0.002375, 0.00475, 0.007125, 0.0095, 0.011875, 0.01425, 0.016625, 0.019, 0.021375, 0.02375, 0.026125, 0.0285, 0.030875, 0.03325])
Ifund_mA = np.array([1.5131, 1.6246, 1.4932, 1.6711, 1.7109, 1.817, 1.7177, 1.6502, 1.6852, 1.7652, 1.8488, 1.9048, 1.9497, 1.9705, 1.9756, 1.9717, 1.9531, 1.9116, 1.8525, 1.7779, 1.6928, 1.6577, 1.7247, 1.8184, 1.7306, 1.6805, 1.5131, 1.6378, 1.5375])
Ifund_A = Ifund_mA/1e3

mid_idx = int(np.floor(len(Idc_A)/2))

# Create array of Pgen and L0 values to try
L0_list = np.linspace(1e-9, 600e-9, 600)
Pgen_list = np.linspace(-6, 0, 13)

t0 = time.time()
ncomplete = 0

# Loop over all L0 values
Iac_results = []
rmse1_results = []
rmse2_results = []
coefs = []
conditions = []

for idxP, Pgen in enumerate(Pgen_list):
	
	for idx, L0 in enumerate(L0_list):
		
		ncomplete += 1
		
		t0_ = time.time()
		
		# Report scan iteration
		logging.main(f"{Fore.LIGHTBLUE_EX}Beginning simulation of L0 = {rd(L0*1e9,1)} nH and Pgen = {rd(Pgen,1)} dBm{standard_color}")
		
		# Prepare simualtion
		lks = LKSystem(Pgen, C_, l_phys, freq, q, L0)
		lks.opt.start_guess_method = GUESS_USE_LAST
		lks.solve(Idc_A, show_plot_on_conv=False)
		
		# Get results
		Iac = np.array([x.Iac for x in lks.solution])
		
		# Calculate error basic error
		rmse1 = np.sqrt( np.mean( np.abs(Ifund_A - Iac)**2 ) )
		
		# Calculate auto-scaled error
		coef = Ifund_A[mid_idx]/Iac[mid_idx]
		Iac_scaled = Iac * coef
		rmse2 = np.sqrt( np.mean( np.abs(Ifund_A - Iac_scaled)**2 ) )
		
		# Add to output data
		Iac_results.append(Iac)
		rmse1_results.append(rmse1)
		rmse2_results.append(rmse2)
		coefs.append(coef)
		conditions.append( (Pgen, L0) )
		
		# Report stats and estimated completion time
		logging.main(f"Standard error assesment: [RMSE = {Fore.LIGHTRED_EX}{rd(rmse1*1e3,3)}m{standard_color}]. Minimum RMSE = {Fore.LIGHTRED_EX}{rd(np.min(rmse1_results)*1e3,3)}m{standard_color}].")
		logging.main(f"Autoscaled error assesment: [RMSE = {Fore.LIGHTRED_EX}{rd(rmse2*1e3,3)}m{standard_color}]. Minimum RMSE = {Fore.LIGHTRED_EX}{rd(np.min(rmse2_results)*1e3,3)}m{standard_color}].")
		
		iter_per_sec = ncomplete/(time.time()-t0)
		t_est = (len(L0_list)*len(Pgen_list)-ncomplete)/iter_per_sec
		logging.main(f"Iteration time = {Fore.LIGHTRED_EX}{rd(time.time()-t0_,1)}{standard_color} sec. Est time remaining = {Fore.LIGHTRED_EX}{rd(t_est, 2)}{standard_color} sec.")

# Find minimum error 
idx_best1 = rmse1_results.index(min(rmse1_results))
idx_best2 = rmse2_results.index(min(rmse2_results))

# Get L0 conditions
L0_best1 = conditions[idx_best1][1]
L0_best2 = conditions[idx_best2][1]

# Get Pg conditions
Pg_best1 = conditions[idx_best1][0]
Pg_best2 = conditions[idx_best2][0]

coef_best = coefs[idx_best2]

plt.plot(Idc_A, Ifund_mA, color='r', label="Measurement", linestyle='dashed', marker='o')
plt.plot(Idc_A, 1e3*Iac_results[idx_best1], color='b', label=f"Sim (L0={rd(L0_best1*1e9,1)} nH) (Pgen = {Pg_best1} dBm)", linestyle='dashed', marker='o')
plt.plot(Idc_A, 1e3*Iac_results[idx_best2]*coef_best, color='g', label=f"Scaled Sim (L0={rd(L0_best2*1e9, 1)} nH) (Pgen = {Pg_best2} dBm)", linestyle='dashed', marker='o')
plt.grid()
plt.xlabel("Bias Current (mA)")
plt.ylabel("AC Current Amplitude (mA)")
plt.title(f"Closest fit")
plt.legend()

print(f"Scaling coefficient used in plot: {coef_best}")

# Save data to Pickle for further analysis
master_data = {"Iac_results":Iac_results, "rmse1_results":rmse1_results, "rmse2_results":rmse2_results, "Idc_A":Idc_A, "coefs":coefs, "conditions":conditions}
with open(filename, 'wb') as f:
	pickle.dump(master_data, f)

plt.show()