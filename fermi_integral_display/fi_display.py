import time
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider

start_time = time.time()

# Constants
kB_eV = 8.617332385*10**-5 # [eV/K]
q = 1.609e-19 #[C]
h = 4.135667e-15 # [eV * s]

#LOAD MATERIAL DATA
data = loadmat('Bi2Te3_QL_scissor_data.mat')
E = np.squeeze(np.array(data['egrid']))
E_min = min(E)
E_max = max(E)
dE = abs(E[1] - E[0])
modes = np.squeeze(np.array(data['M']))
VELO = np.squeeze(np.array(data['VELO']))
ones = np.ones_like(E)
thickness = np.squeeze(np.array(data['thickness']))

#Transmission functions to compare
transmission_MFP = modes       #Constant mean free path
transmission_TAU = modes*VELO  #Constant relaxation time
transmission_DOS = VELO*VELO   #DOS scattering

mean_mfp = 25e-9
T_min = 200
T_max = 2000
T_step = 100
T_vec = np.arange(T_min,T_max, T_step)

#Scaling Factors 
f1_scale = 0.5
f2_scale = 0.5
 
moment0_MFP = np.zeros([len(E), len(T_vec)])
moment1_MFP = np.zeros([len(E), len(T_vec)])
moment2_MFP = np.zeros([len(E), len(T_vec)])
moment0_TAU = np.zeros([len(E), len(T_vec)])
moment1_TAU = np.zeros([len(E), len(T_vec)])
moment2_TAU = np.zeros([len(E), len(T_vec)])
moment0_DOS = np.zeros([len(E), len(T_vec)])
moment1_DOS = np.zeros([len(E), len(T_vec)])
moment2_DOS = np.zeros([len(E), len(T_vec)])
tau_values = np.zeros([len(E), len(T_vec)])
K_elph_values = np.zeros([len(E), len(T_vec)])   

#FUNCTIONS---------------------------------------------------------------------------------------------------------------------------------------------------------
def fermi(E,mu,T):
    #Fermi Window (Zeroth Order)
    return (1/(kB_eV*T))*np.exp((E-mu*ones)/(kB_eV*T))/(ones +  np.exp((E-mu*ones)/(kB_eV*T)) )**2

def fermi1(E,mu,T):
    #First Order Fermi Window
    return fermi(E,mu, T)*(E-mu*ones)/(kB_eV*T) 

def fermi2(E,mu,T):
    #Second Order Fermi Window
    return fermi(E,mu,T)*((E-mu*ones)/(kB_eV*T))**2

def constant_calc(mu, temp, transmission, avg_mfp):
    "Returns constant that fixes conductance between different scattering models"
    f0 = fermi(E, mu, temp)
    numerator = np.trapz(transmission*f0, E)
    denominator = np.trapz(modes*f0, E)
    ratio = numerator/denominator 
    return avg_mfp/ratio

#CALCULATION------------------------------------------------------------------------------------------------------------------------------------------------------
print('CALCULATING THERMOELECTRIC MOMENTS')
for k in range(len(E)):
    percent = 100*k/len(E)
    print('Completed: ', "%.2f" % percent, end="\r")
    for j in range(len(T_vec)):
        f0 = fermi(E, E[k], T_vec[j])
        f1 = fermi1(E, E[k], T_vec[j])
        f2 = fermi2(E, E[k], T_vec[j])
        tau_values[k, j] = constant_calc(E[k], T_vec[j], transmission_TAU, mean_mfp)
        K_elph_values[k, j] = constant_calc(E[k], T_vec[j], transmission_DOS, mean_mfp)

        moment0_MFP[k, j] = np.trapz(f0*transmission_MFP,E)
        moment1_MFP[k, j] = np.trapz(f1*transmission_MFP, E)
        moment2_MFP[k, j] = np.trapz(f2*transmission_MFP, E)

        moment0_TAU[k, j] = np.trapz(f0*transmission_TAU,E)
        moment1_TAU[k, j] = np.trapz(f1*transmission_TAU, E)
        moment2_TAU[k, j] = np.trapz(f2*transmission_TAU, E)
       
        moment0_DOS[k, j] = np.trapz(f0*transmission_DOS,E)
        moment1_DOS[k, j] = np.trapz(f1*transmission_DOS, E)
        moment2_DOS[k, j] = np.trapz(f2*transmission_DOS, E)

print(' ', end="/r")
print('100 % complete')
print('Calculated thermoelectric moments for ', len(T_vec), ' temperatures on an energy grid with a resolution of ', dE, ' eV.')

fig1, axs = plt.subplots(3,3)

#INITIALIZATION--------------------------------------------------------------------------------------------------------------------------------------------------------------
initial_T = 300
initial_T_index = int((initial_T-T_min)/T_step)
initial_mu = 0
initial_mu_index = int((initial_mu-E_min)/dE)

f0 = fermi(E, initial_mu, initial_T)
f1 = fermi1(E, initial_mu, initial_T)
f2 = fermi2(E, initial_mu, initial_T)

#Figure [0,0]
T0_MFP, = axs[0,0].plot(E,mean_mfp*transmission_MFP,'b', lw=2)
T0_TAU, = axs[0,0].plot(E, tau_values[initial_mu_index, initial_T_index]*transmission_TAU, 'r', lw=2)
T0_DOS, = axs[0,0].plot(E, K_elph_values[initial_mu_index, initial_T_index]*transmission_DOS, 'k', lw=2)
F0, = axs[0,0].plot(E, f0,'g--', lw=2)
axs[0,0].legend(['Constant MFP', 'Constant Relaxation Time', 'DOS Scattering', 'Fermi Window'])
#Figure [1,0]
F1, = axs[1,0].plot(E, f1_scale*f1,'g--', lw=2)
T1_MFP, = axs[1,0].plot(E,mean_mfp*transmission_MFP,'b', lw=2)
T1_TAU, = axs[1,0].plot(E, tau_values[initial_mu_index, initial_T_index]*transmission_TAU,'r', lw=2)
T1_DOS, = axs[1,0].plot(E, K_elph_values[initial_mu_index, initial_T_index]*transmission_DOS,'k', lw=2)
#Figure [2,0]
F2, = axs[2,0].plot(E, f2_scale*f2,'g--', lw=2)
T2_MFP, = axs[2,0].plot(E,mean_mfp*transmission_MFP,'b', lw=2)
T2_TAU, = axs[2,0].plot(E, tau_values[initial_mu_index, initial_T_index]*transmission_TAU,'r', lw=2)
T2_DOS, = axs[2,0].plot(E, K_elph_values[initial_mu_index, initial_T_index]*transmission_DOS, 'k',lw=2)

#Figure [0,1]
I0_MFP, = axs[0,1].plot(E, mean_mfp*transmission_MFP*f0, 'b')
I0_TAU, = axs[0,1].plot(E, tau_values[initial_mu_index, initial_T_index]*transmission_TAU*f0, 'r')
I0_DOS, = axs[0,1].plot(E, K_elph_values[initial_mu_index, initial_T_index]*transmission_DOS*f0, 'k')
#Figure [1,1]
I1_MFP, = axs[1,1].plot(E, mean_mfp*transmission_MFP*f1, 'b')
I1_TAU, = axs[1,1].plot(E, tau_values[initial_mu_index, initial_T_index]*transmission_TAU*f1, 'r')
I1_DOS, = axs[1,1].plot(E, K_elph_values[initial_mu_index, initial_T_index]*transmission_DOS*f1, 'k')
#Figure [2,1]
I2_MFP, = axs[2,1].plot(E, mean_mfp*transmission_MFP*f2, 'b')
I2_TAU, = axs[2,1].plot(E, tau_values[initial_mu_index, initial_T_index]*transmission_TAU*f2, 'r')
I2_DOS, = axs[2,1].plot(E, K_elph_values[initial_mu_index, initial_T_index]*transmission_DOS*f2, 'k')

#Figure [0,2]
M0_MFP, = axs[0,2].plot(E, moment0_MFP[:, initial_T_index], 'b')
M0_TAU, = axs[0,2].plot(E, moment0_TAU[:, initial_T_index], 'r')
M0_DOS, = axs[0,2].plot(E, moment0_DOS[:, initial_T_index], 'k')

V0, = axs[0,2].plot([initial_mu, initial_mu], [0,50], '--')
H0_MFP, = axs[0,2].plot([-1.5, 1.5], [max(moment0_MFP[:, initial_T_index]), max(moment0_MFP[:, initial_T_index])], 'b--')
H0_TAU, = axs[0,2].plot([-1.5, 1.5], [max(moment0_TAU[:, initial_T_index]), max(moment0_TAU[:, initial_T_index])], 'r--')
H0_DOS, = axs[0,2].plot([-1.5, 1.5], [max(moment0_DOS[:, initial_T_index]), max(moment0_DOS[:, initial_T_index])], 'k--')

#Figure [1,2]
M1_MFP, = axs[1,2].plot(E, moment1_MFP[:, initial_T_index ], 'b') 
M1_TAU, = axs[1,2].plot(E, moment1_TAU[:, initial_T_index], 'r')
M1_DOS, = axs[1,2].plot(E, moment1_DOS[:, initial_T_index], 'k')

V1, = axs[1,2].plot([initial_mu, initial_mu], [0,50], '--')
H1_MFP, = axs[1,2].plot([-1.5, 1.5], [max(moment1_MFP[:, initial_T_index]), max(moment1_MFP[:, initial_T_index])], 'b--')
H1_TAU, = axs[1,2].plot([-1.5, 1.5], [max(moment1_TAU[:, initial_T_index]), max(moment1_TAU[:, initial_T_index])], 'r--')
H1_DOS, = axs[1,2].plot([-1.5, 1.5], [max(moment1_DOS[:, initial_T_index]), max(moment1_DOS[:, initial_T_index])], 'k--')

#Figure [2,2]
M2_MFP, = axs[2,2].plot(E, moment2_MFP[:, initial_T_index ], 'b') 
M2_TAU, = axs[2,2].plot(E, moment2_TAU[:, initial_T_index], 'r')
M2_DOS, = axs[2,2].plot(E, moment2_DOS[:, initial_T_index], 'k')

V2, = axs[2,2].plot([initial_mu, initial_mu], [0,50], '--')
H2_MFP, = axs[2,2].plot([-1.5, 1.5], [max(moment2_MFP[:, initial_T_index]), max(moment2_MFP[:, initial_T_index])], 'b--')
H2_TAU, = axs[2,2].plot([-1.5, 1.5], [max(moment2_TAU[:, initial_T_index]), max(moment2_TAU[:, initial_T_index])], 'r--')
H2_DOS, = axs[2,2].plot([-1.5, 1.5], [max(moment2_DOS[:, initial_T_index]), max(moment2_DOS[:, initial_T_index])], 'k--')

#SLIDERS------------------------------------------------------------------------------------------------------------------------------------------------------------------
mu_loc= plt.axes([0.25, 0.025, 0.50, 0.02])
T_loc = plt.axes([0.25, 0.00, 0.50, 0.02])
mu_amp = Slider(mu_loc, 'Chemical Potential', -1.5, 1.5, valinit=initial_mu)
T_amp = Slider(T_loc, 'Temperature', T_min, T_max, valinit=initial_T)

#LABELS---------------------------------------------------------------------------------------------------------------------------------------------------------------------

#First Column
axs[2,0].set_xlabel('Energy Level (eV)')
axs[0,0].set_title('Zeroth Order')
axs[1,0].set_title('First Order')
axs[2,0].set_title('Second Order')

#Second Column
axs[2,1].set_xlabel('Energy Level (eV)')
axs[0,1].set_title('I_0 Integrand')
axs[1,1].set_title('I_1 Integrand')
axs[2,1].set_title('I_2 Integrand')
axs[0,1].set_xlim([-1.0, 1.0])
axs[1,1].set_xlim([-1.0, 1.0])
axs[2,1].set_xlim([-1.0, 1.0])

#Third Column
axs[2,2].set_xlabel('Fermi Level (eV)')
axs[0,2].set_ylabel('G [S/m]')
axs[0,2].set_title('Conductance')
axs[1,2].set_title('Seebeck Coefficient')
axs[2,2].set_title('Thermal Conductivity')
axs[1,2].set_ylabel('S [\mu V/K]')
axs[2,2].set_ylabel('\kappa_e [W/(m*K)]')
axs[0,2].set_xlim([-1.0, 1.0])
axs[1,2].set_xlim([-1.0, 1.0])
axs[2,2].set_xlim([-1.0, 1.0])

def update(val): 
    mu = mu_amp.val
    T = T_amp.val
    mu_index = int((mu-E_min)/dE)
    mu = E[mu_index]
    T_index = int((T-T_min)/T_step)
    T = T_vec[T_index]
  
    f0 = fermi(E, mu, T)
    f1 = fermi1(E, mu, T)
    f2 = fermi2(E, mu, T)
    m0_MFP = transmission_MFP*f0
    m1_MFP = transmission_MFP*f1
    m2_MFP = transmission_MFP*f2
    m0_TAU = transmission_TAU*f0
    m1_TAU = transmission_TAU*f1
    m2_TAU = transmission_TAU*f2
    m0_DOS = transmission_DOS*f0
    m1_DOS = transmission_DOS*f1
    m2_DOS = transmission_DOS*f2
#FIRST COLUMN
    F0.set_ydata(f0*0.1*max(max(mean_mfp*transmission_MFP), max(tau_values[mu_index, T_index]*transmission_TAU), max(K_elph_values[mu_index, T_index]*transmission_DOS)))
    F1.set_ydata(f1_scale*f1*0.1*max(max(mean_mfp*transmission_MFP), max(tau_values[mu_index, T_index]*transmission_TAU), max(K_elph_values[mu_index, T_index]*transmission_DOS)))
    F2.set_ydata(f2_scale*f2*0.1*max(max(mean_mfp*transmission_MFP), max(tau_values[mu_index, T_index]*transmission_TAU), max(K_elph_values[mu_index, T_index]*transmission_DOS)))
    axs[0,0].set_ylim(0.0, 1.1*max(max(mean_mfp*transmission_MFP), max(tau_values[mu_index, T_index]*transmission_TAU), max(K_elph_values[mu_index, T_index]*transmission_DOS)))
    axs[1,0].set_ylim(-f1_scale*0.8*max(max(mean_mfp*transmission_MFP), max(tau_values[mu_index, T_index]*transmission_TAU), max(K_elph_values[mu_index, T_index]*transmission_DOS)), 1.1*max(max(mean_mfp*transmission_MFP), max(tau_values[mu_index, T_index]*transmission_TAU), max(K_elph_values[mu_index, T_index]*transmission_DOS)))
    axs[2,0].set_ylim(0.0, max(max(mean_mfp*transmission_MFP),max(tau_values[mu_index, T_index]*transmission_TAU), max(K_elph_values[mu_index, T_index]*transmission_DOS)))
    T0_MFP.set_ydata(mean_mfp*transmission_MFP)
    T1_MFP.set_ydata(mean_mfp*transmission_MFP)
    T2_MFP.set_ydata(mean_mfp*transmission_MFP)
    T0_TAU.set_ydata(tau_values[mu_index, T_index]*transmission_TAU)
    T1_TAU.set_ydata(tau_values[mu_index, T_index]*transmission_TAU)
    T2_TAU.set_ydata(tau_values[mu_index, T_index]*transmission_TAU)
    T0_DOS.set_ydata(K_elph_values[mu_index, T_index]*transmission_DOS)
    T1_DOS.set_ydata(K_elph_values[mu_index, T_index]*transmission_DOS)
    T2_DOS.set_ydata(K_elph_values[mu_index, T_index]*transmission_DOS)
#SECOND COLUMN
    I0_MFP.set_ydata(mean_mfp*m0_MFP)
    I0_TAU.set_ydata(tau_values[mu_index, T_index]*m0_TAU)
    I0_DOS.set_ydata(K_elph_values[mu_index, T_index]*m0_DOS)
    axs[0,1].set_ylim([0,1.1*max(max(mean_mfp*m0_MFP), max(tau_values[mu_index, T_index]*m0_TAU), max(K_elph_values[mu_index, T_index]*m0_DOS))])
    I1_MFP.set_ydata(mean_mfp*m1_MFP)
    I1_TAU.set_ydata(tau_values[mu_index, T_index]*m1_TAU)
    I1_DOS.set_ydata(K_elph_values[mu_index, T_index]*m1_DOS)
    axs[1,1].set_ylim([1.1*min(min(mean_mfp*m1_MFP), min(tau_values[mu_index, T_index]*m1_TAU), min(K_elph_values[mu_index, T_index]*m1_DOS)), 1.1*max(max(mean_mfp*m1_MFP), max(tau_values[mu_index, T_index]*m1_TAU), max(K_elph_values[mu_index, T_index]*m1_DOS))])
    I2_MFP.set_ydata(mean_mfp*m2_MFP)
    I2_TAU.set_ydata(tau_values[mu_index, T_index]*m2_TAU)
    I2_DOS.set_ydata(K_elph_values[mu_index, T_index]*m2_DOS)
    axs[2,1].set_ylim([0, 1.1*max(max(mean_mfp*m2_MFP), max(tau_values[mu_index,T_index]*m2_TAU), max(K_elph_values[mu_index, T_index]*m2_DOS))])
#THIRD COLUMN 
    m0_MFP = moment0_MFP[:,T_index] 
    G_MFP = 2*q*q*m0_MFP/h
    G_MFP = 25e-9*G_MFP/(q*thickness)
    M0_MFP.set_ydata(G_MFP)
    G_MFP_mu = G_MFP[mu_index]
    H0_MFP.set_data([-1.5, 1.5], [G_MFP_mu, G_MFP_mu])
    
    m0_TAU = moment0_TAU[:, T_index]
    G_TAU = 2*q*q*m0_TAU/h
    G_TAU = tau_values[:, T_index]*G_TAU/(q*thickness)
    M0_TAU.set_ydata(G_TAU)
    G_TAU_mu = G_TAU[mu_index]
    H0_TAU.set_data([-1.5, 1.5], [G_TAU_mu, G_TAU_mu])

    m0_DOS = moment0_DOS[:, T_index]
    G_DOS = 2*q*q*m0_DOS/h
    G_DOS = K_elph_values[:, T_index]*G_DOS/(q*thickness)
    M0_DOS.set_ydata(G_DOS)
    G_DOS_mu = G_DOS[mu_index]
    H0_DOS.set_data([-1.5, 1.5], [G_DOS_mu, G_DOS_mu])

    axs[0,2].set_ylim([0, 1.1*max(max(G_MFP), max(G_TAU), max(G_DOS))])
    V0.set_data([mu, mu], [0, 1.2*max(max(G_MFP), max(G_TAU), max(G_DOS))])

    m1_MFP = moment1_MFP[:, T_index]
    S_MFP = -kB_eV*m1_MFP/m0_MFP
    S_MFP = 1e6*S_MFP
    S_MFP_mu = S_MFP[mu_index]
    H1_MFP.set_data([-1.5, 1.5], [S_MFP_mu, S_MFP_mu])
    M1_MFP.set_ydata(S_MFP)
   
    m1_TAU = moment1_TAU[:, T_index]
    S_TAU = -kB_eV*m1_TAU/m0_TAU
    S_TAU = 1e6*S_TAU
    S_TAU_mu = S_TAU[mu_index]
    H1_TAU.set_data([-1.5, 1.5], [S_TAU_mu, S_TAU_mu])
    M1_TAU.set_ydata(S_TAU)

    m1_DOS = moment1_DOS[:, T_index]
    S_DOS = -kB_eV*m1_DOS/m0_DOS
    S_DOS = 1e6*S_DOS
    S_DOS_mu = S_DOS[mu_index]
    H1_DOS.set_data([-1.5, 1.5], [S_DOS_mu, S_DOS_mu])
    M1_DOS.set_ydata(S_DOS)

    axs[1,2].set_ylim([1.1*min(min(S_MFP), min(S_TAU), min(S_DOS)), 1.1*max(max(S_MFP), max(S_TAU), max(S_DOS))])
    V1.set_data([mu, mu], [1.2*min(min(S_MFP), min(S_TAU), min(S_DOS)), 1.2*max(max(S_MFP), max(S_TAU), max(S_DOS))])

    m2_MFP = moment2_MFP[:, T_index]
    k_e_MFP = (2*T*q*kB_eV*kB_eV/h)*(m2_MFP - m1_MFP*m1_MFP/m0_MFP)
    k_e_MFP = mean_mfp*k_e_MFP
    k_e_MFP = 1.5*ones + k_e_MFP/thickness
    k_e_MFP_mu = k_e_MFP[mu_index]
    H2_MFP.set_data([-1.5, 1.5], [k_e_MFP_mu, k_e_MFP_mu])
    M2_MFP.set_ydata(k_e_MFP)

    m2_TAU = moment2_TAU[:, T_index]
    k_e_TAU = (2*T*q*kB_eV*kB_eV/h)*(m2_TAU - m1_TAU*m1_TAU/m0_TAU)
    k_e_TAU = tau_values[:, T_index]*k_e_TAU
    k_e_TAU = 1.5*ones + k_e_TAU/thickness
    k_e_TAU_mu = k_e_TAU[mu_index]
    H2_TAU.set_data([-1.5, 1.5], [k_e_TAU_mu, k_e_TAU_mu])
    M2_TAU.set_ydata(k_e_TAU)
   
    m2_DOS = moment2_DOS[:, T_index]
    k_e_DOS = (2*T*q*kB_eV*kB_eV/h)*(m2_DOS - m1_DOS*m1_DOS/m0_DOS)
    k_e_DOS = K_elph_values[:, T_index]*k_e_DOS
    k_e_DOS = 1.5*ones + k_e_DOS/thickness
    k_e_DOS_mu = k_e_DOS[mu_index]
    H2_DOS.set_data([-1.5, 1.5], [k_e_DOS_mu, k_e_DOS_mu])
    M2_DOS.set_ydata(k_e_DOS)
    
    axs[2,2].set_ylim([0, 1.1*max(max(k_e_MFP), max(k_e_TAU), max(k_e_DOS))])
    V2.set_data([mu, mu], [0, 1.2*max(max(k_e_MFP), max(k_e_TAU), max(k_e_DOS))])

# redraw canvas while idle
    fig1.canvas.draw_idle()

elapsed_time = time.time() - start_time
print('Time taken: ', elapsed_time, 'seconds')
mu_amp.on_changed(update)
T_amp.on_changed(update)
plt.show()





