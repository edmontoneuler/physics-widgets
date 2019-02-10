import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
from matplotlib.widgets import Slider


#k grid
theta_min = -np.pi
theta_max = np.pi
theta = np.array(np.linspace(theta_min, theta_max, 1000))


def free_energy(theta, B, kappa=1):
    return -1*B*np.cos(theta) -kappa*np.cos(theta)**2

B_init = 0

fig, axs = plt.subplots()
f, = axs.plot(np.cos(theta), free_energy(theta, B_init), 'k', lw=2)

B_loc = plt.axes([0.25, 0.01, 0.50, 0.02])
B_amp = Slider(B_loc, 'B (T) ', -3, 3, valinit = B_init)

axs.set_xlabel('Magnetization Angle (cos(rad))')
axs.set_ylabel('Free Energy')
axs.set_title('Hysteresis Model', fontsize=20)

#axs.set_ylim([-0.25, 0.75])

def update(val):
    B = B_amp.val
    F = free_energy(theta,B)
    f.set_ydata(F)
    axs.set_ylim([min(F), max(F) +0.5])

fig.canvas.draw_idle()
B_amp.on_changed(update)
plt.show()
