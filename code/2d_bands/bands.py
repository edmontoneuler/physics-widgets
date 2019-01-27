import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

x = np.linspace(-3, 3, 20)
y = x 
X, Y = np.meshgrid(x,y)

def three_param_band(A, B, C, gap=0, xgrid = X, ygrid = Y):
    return A*X*X + B*X*Y + C*Y*Y + gap*np.ones_like(X)

val = three_param_band(-0.2, 0, -0.5)
cond = three_param_band(5, 0, 5, gap=1)

V, = ax.plot_surface(X, Y, val)
C, = ax.plot_surface(X, Y, cond)

A_loc = plt.axes([0.25, -0.00, 0.50, 0.02])
A_amp = Slider(A_loc, 'A', -1, 1, valinit = 
plt.show()
