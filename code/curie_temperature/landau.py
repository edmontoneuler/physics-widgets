import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
from matplotlib.widgets import Slider


#M grid
M_min = -1
M_max = 1
M = np.linspace(M_min, M_max, 1000)


def free_energy(M, alpha, beta):
    return alpha*M**2 + beta*M**4

alpha_init = 2
beta_init = 1

fig, axs = plt.subplots()
f, = axs.plot(M, free_energy(M, alpha_init, beta_init), 'k', lw=2)

alpha_loc = plt.axes([0.25, 0.01, 0.50, 0.02])
alpha_amp = Slider(alpha_loc, 'A ', -2, 2, valinit = alpha_init)

beta_loc = plt.axes([0.25, 0.04, 0.50, 0.02])
beta_amp = Slider(beta_loc,'beta', 0,10 , valinit = beta_init)

axs.set_xlabel('M/M_s')
axs.set_ylabel('Free Energy (eV)')
axs.set_title('Paramagnetic-Ferromagnetic Transition', fontsize=20)

#axs.set_ylim([-0.25, 0.75])

def update(val):
    alpha = alpha_amp.val
    beta = beta_amp.val
    F = free_energy(M, alpha, beta)
    f.set_ydata(F)
    axs.set_ylim([min(F) - 0.5, max(F) +0.5])

fig.canvas.draw_idle()
alpha_amp.on_changed(update)
beta_amp.on_changed(update)
plt.show()
