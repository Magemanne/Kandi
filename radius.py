#original found in https://scicomp.stackexchange.com/questions/34084/numerical-solution-to-rayleigh-plesset-equation-in-python
#changed the equation which is used.
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint

# parameters

rho = 1000
sigma = 0.07169
miu = 8.9*10**(-4)
kv = 1.206*10**(-7)
P_g = 0
P_0 = 100000
k = 1.33
fs = 10**(-15)
A = 10**(-10)

# define equations
def equation(y0, t):
    R, u = y0
    #return u, (P_g-P_0-2*sigma/R-4*miu*u/R+(2*sigma/R_0+P_0-P_g)*(R_0/R)**(3*k))/(R*rho)-3*u**2/(2*R)
    return u,(((P_g-P_0)/rho - 3/2*u**2 - 4*kv*u/R - 2*sigma/(R*rho))/R)
# initial conditions
R_0 = 34.424*A
u_0 = 0



time = np.arange(0, 19100, 1)
time = time*fs

R_1 = odeint(equation, [R_0, u_0], time)

V = R_1[:,1]
R = R_1[:,0]
mtimes = time
K = np.array((mtimes,V,R)).T
np.savetxt('test.out', K, delimiter=" ",newline='\n',header="data" )
#plot results

fig, ax1 = plt.subplots()

ax1.set_xlabel("T/$\mu$s")
ax1.set_ylabel("R/$\mu$m", color = "red")
ax1.plot(mtimes, R, linewidth = 0.7, label = "Bubble Radius", color =
"red")

ax2 = ax1.twinx()

ax2.set_ylabel("$\dot{R}$/$ms^{-1}$")
ax2.plot(mtimes, V, linestyle = "dashed", color = "black", linewidth =
0.7, label = "Radial Velocity Bubble")
ax1.legend()
ax2.legend(loc = "lower right")

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
