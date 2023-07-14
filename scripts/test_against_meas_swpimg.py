from Lkpy.core import *
import time
import pickle
import os

Pgen = 0 # dBm
C_ = 150.2e-12
l_phys = 0.5
freq = 10e9
# L0 = 850e-9*l_phys
# qs = [0.185, 0.19, 0.195, 0.20, 0.205, 0.21]
q = 0.19
L0s = [750e-9, 800e-9, 850e-9, 900e-9, 950e-9, 1000e-9]

filename = "last_simulation.pkl"

# Open S21 data File
with open("cryostat_sparams.pkl", 'rb') as fh:
	S21_data = pickle.load(fh)

# Get target values from MATLAB
Idc_A = np.array([-0.03325, -0.030875, -0.0285, -0.026125, -0.02375, -0.021375, -0.019, -0.016625, -0.01425, -0.011875, -0.0095, -0.007125, -0.00475, -0.002375, 0, 0.002375, 0.00475, 0.007125, 0.0095, 0.011875, 0.01425, 0.016625, 0.019, 0.021375, 0.02375, 0.026125, 0.0285, 0.030875, 0.03325])
Ifund_mA = np.array([1.5131, 1.6246, 1.4932, 1.6711, 1.7109, 1.817, 1.7177, 1.6502, 1.6852, 1.7652, 1.8488, 1.9048, 1.9497, 1.9705, 1.9756, 1.9717, 1.9531, 1.9116, 1.8525, 1.7779, 1.6928, 1.6577, 1.7247, 1.8184, 1.7306, 1.6805, 1.5131, 1.6378, 1.5375])
Ifund_A = Ifund_mA/1e3

mid_idx = int(np.floor(len(Idc_A)/2))

# for q in qs:
for L0 in L0s:
	
	logging.info(f"Beginning sweep for q={Fore.GREEN}{rd(q*1e3)}{Style.RESET_ALL} mA")
	
	# Initialize system
	lks = LKSystem(Pgen, C_, l_phys, freq, q, L0)

	# Change options (applies to all simulators)
	lks.setopt('start_guess_method', GUESS_USE_LAST)
	lks.setopt('max_iter', 200)
	lks.configure_loss(sparam_data=S21_data)

	# Select simulator for system to use
	lks.select_simulator(SIMULATOR_ABCD)

	# Run hybrid simulation
	lks.solve(Idc_A)
	lks.solve(Idc_A, simulator=SIMULATOR_P0)

	IL_mA = lks.get_solution(simulator=SIMULATOR_ABCD, parameter='IL_w', element=0)*1e3
	IL_mAp0 = lks.get_solution(simulator=SIMULATOR_P0, parameter='Ig_w', element=0)*1e3
	Idc_ABCD = lks.get_solution(simulator=SIMULATOR_ABCD, parameter='Ibias_c')*1e3
	Idc_P0 = lks.get_solution(simulator=SIMULATOR_P0, parameter='Ibias_c')*1e3
	
	plt.cla()
	plt.plot(Idc_A*1e3, Ifund_mA, color=(0.7, 0, 0), label="Measurement", linestyle='dashed', marker='o')
	plt.plot(Idc_ABCD, IL_mA, color=(0, .7, 0.4), label="ABCD Simulation", linestyle='dashed', marker='o')
	plt.plot(Idc_P0, IL_mAp0, color=(0, 0.4, 0.7), label="P0 Simulation", linestyle='dashed', marker='o')
	plt.grid()
	plt.legend()
	plt.xlabel("Bias Current (mA)")
	plt.ylabel("AC Current Amplitude (mA)")
	# plt.title(f"Fundamental at VNA, q={rd(q*1e3)}mA")
	plt.title(f"Fundamental at VNA, L0={rd(L0*1e9)} nH")
	plt.legend()
	
	# Force0y
	ax = plt.gca()
	yl = ax.get_ylim()
	ax.set_ylim([0, yl[1]])
	
	# plt.savefig(os.path.join("q sweep", f"plot_q={rd(q*1e3)}mA.png"))
	plt.savefig(os.path.join("L0 sweep", f"plot_L0={rd(L0*1e9)}nH.png"))


