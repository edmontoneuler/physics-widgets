import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
from matplotlib.widgets import Slider

#The three-band Kane model, derived via k.p perturbation theory, is a 
#simple model that describes the bandstructure (for small k) of the lower most conduction band and two uppermost valence bands in a semiconductor, assuming that the energies and wave functions of these electrons are known at a single point in reciprocal space. 

#This applet plots the three band structures of this model as a function of band gap
#and P parameter, which can be controlled by sliders to provide slightly more intuition for their physical effect

#k grid
a = 1e-9 # lattice parameter [m]
k_min = -2*np.pi/a
k_max = 2*np.pi/a
k_vec = np.linspace(k_min, k_max, 1000)
ones = np.ones_like(k_vec)

hbar = 1.0545718e-34 # [J * s]
m0 = 9.109384e-31 # [kg]
q = 1.602e-19     # [C]

def valence_one(k_vec):
    return -hbar*hbar*k_vec*k_vec/(2*m0*q)

def valence_two(k_vec, gap, P):
    ones = np.ones_like(k_vec)
    return 0.5*gap*ones - hbar*hbar*k_vec*k_vec/(2*m0*q) - 0.5*np.sqrt(gap*gap*ones +8*hbar*hbar*k_vec*k_vec*P*P/(m0*m0))

def conduction_one(k_vec, gap, P):
    ones = np.ones_like(k_vec)
    return 0.5*gap*ones + hbar*hbar*k_vec*k_vec/(2*m0*q) + 0.5*np.sqrt(gap*gap*ones +8*hbar*hbar*k_vec*k_vec*P*P/(m0*m0))

#Initialization
fig, axs = plt.subplots()

P_init = 1e-8
gap_init = 0.5

c1, = axs.plot(k_vec, conduction_one(k_vec, gap_init, P_init), 'r', lw=2)
v1, = axs.plot(k_vec, valence_one(k_vec), 'b', lw=2)
v2, = axs.plot(k_vec, valence_two(k_vec, gap_init, P_init), 'g', lw=2)

#Sliders
P_loc = plt.axes([0.25, 0.00, 0.50, 0.02])
P_amp = Slider(P_loc, 'P', 0, 1e-5, valinit = P_init)

gap_loc = plt.axes([0.25, 0.05, 0.50, 0.02])
gap_amp = Slider(gap_loc,'Band gap (eV)', 0.05, 1.0, valinit = gap_init)

#Labels
axs.set_xlabel('k (m^-1)')
axs.set_ylabel('E (eV)')
axs.legend(['C1', 'V1', 'V2'])

def update(val):
    P = P_amp.val
    P = P
    gap = gap_amp.val
    gap = gap
    
    cond1 = conduction_one(k_vec, gap, P)
    val1 = valence_one(k_vec)
    val2 = valence_two(k_vec, gap, P)
    
    c1.set_ydata(cond1)
    v1.set_ydata(val1)
    v2.set_ydata(val2)


fig.canvas.draw_idle()
P_amp.on_changed(update)
gap_amp.on_changed(update)

plt.show()
