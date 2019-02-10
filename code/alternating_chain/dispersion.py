import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
from matplotlib.widgets import Slider


#k grid
a = 1e-9
k_min = -np.pi/a
k_max = np.pi/a
kgrid = np.array(np.linspace(k_min, k_max, 1000))


def acoustic(kgrid, k1, k2, a = 1e-9):
    ones = np.ones_like(kgrid)
    radical = np.sqrt(k1*k1*ones +k2*k2*ones + 2*k1*k2*np.cos(a*kgrid))
    return np.sqrt( (k1+k2)*ones - radical )
def optical(kgrid, k1, k2, a=1e-9):
    ones = np.ones_like(kgrid)
    radical = np.sqrt(k1*k1*ones +k2*k2*ones + 2*k1*k2*np.cos(a*kgrid))
    return np.sqrt( (k1+k2)*ones + radical )

k1_init = 1
k2_init = 2

fig, axs = plt.subplots()
a, = axs.plot(kgrid, acoustic(kgrid, k1_init, k2_init), 'k', lw=2)
o, = axs.plot(kgrid, optical(kgrid, k1_init, k2_init), 'k', lw=2)

k1_loc = plt.axes([0.25, 0.01, 0.50, 0.02])
k1_amp = Slider(k1_loc, 'k1 ', 0.2, 1, valinit = k1_init)

k2_loc = plt.axes([0.25, 0.04, 0.50, 0.02])
k2_amp = Slider(k2_loc,'k2', 0.2,1 , valinit = k2_init)

axs.set_xlabel('k')
axs.set_ylabel('Omega')
axs.set_title('Diatomic Chain Dispersion', fontsize=20)

#axs.set_ylim([-0.25, 0.75])

def update(val):
    k1 = k1_amp.val
    k2 = k2_amp.val
    O = optical(kgrid, k1, k2)
    A = acoustic(kgrid, k1, k2)
    o.set_ydata(O)
    a.set_ydata(A)
    axs.set_ylim([0, max(O) +0.5])

fig.canvas.draw_idle()
k1_amp.on_changed(update)
k2_amp.on_changed(update)
plt.show()
