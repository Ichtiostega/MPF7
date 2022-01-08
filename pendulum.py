from math import pi, sin
from numpy.lib.function_base import average
from scipy.integrate import odeint
import numpy as np
from tabulate import tabulate
from matplotlib import pyplot as plt

dt = 0.1
g = 10
l = 10
steps = 60

class Euler:
    def __init__(self, omega, theta):
        self.omega = omega
        self.theta = theta

    
    def iterate(self):
        omega = self.omega - (g/l)*sin(self.theta)*dt
        theta = self.theta + self.omega * dt
        self.omega, self.theta = omega, theta

    def iterate_nosin(self):
        omega = self.omega - (g/l)*self.theta*dt
        theta = self.theta + self.omega * dt
        self.omega, self.theta = omega, theta

    def __str__(self) :
        return f'omega:\t{self.omega}\ttheta:\t{self.theta}'


class EulerCromer:
    def __init__(self, omega, theta):
        self.omega = omega
        self.theta = theta

    
    def iterate(self):
        omega = self.omega - (g/l)*sin(self.theta)*dt
        theta = self.theta + omega * dt
        self.omega, self.theta = omega, theta

    def iterate_nosin(self):
        omega = self.omega - (g/l)*self.theta*dt
        theta = self.theta + omega * dt
        self.omega, self.theta = omega, theta

    def __str__(self) :
        return f'omega:\t{self.omega}\ttheta:\t{self.theta}'


class ODE:
    def __init__(self, theta, omega):
        self.omega = omega
        self.theta = theta

    @staticmethod
    def _diff(y, t):
        theta, omega = y
        return omega, -(g/l)*sin(theta)

    @staticmethod
    def _diff_nosin(y, t):
        theta, omega = y
        return omega, -(g/l)*theta

    def get_sim(self, iter):
        out = odeint(ODE._diff, [self.theta, self.omega], [i*dt for i in range(0, iter)])
        return out

    def get_sim_nosin(self, iter):
        out = odeint(ODE._diff_nosin, [self.theta, self.omega], [i*dt for i in range(0, iter)])
        return out

theta = 2*pi/10000
ode = ODE(0, theta)
ode_out = ode.get_sim(steps)

data = []
e = Euler(0, theta)
ec = EulerCromer(0, theta)
for i in range(steps):
    data.append([e.omega, e.theta, ec.omega, ec.theta, -ode_out[i][0], ode_out[i][1]])
    e.iterate()
    ec.iterate()

# print(tabulate(data, headers=['e_o', 'e_t', 'ec_o', 'ec_t', 'ode_o', 'ode_t']))

x = [i*dt for i in range(0, steps)]

thetas = list(zip(*data))[1::2]
omegas = list(zip(*data))[::2]

plt.plot(x[::3], thetas[0][::3], '.', x[1::3], thetas[1][1::3], '.', x[2::3], thetas[2][2::3], '.')
pos = plt.gca().get_position()
pos.x0 = 0.15
plt.gca().set_position(pos)
plt.title(f'Wartości odchylenia dla metod')
plt.xlabel('t[s]')
plt.ylabel('Theta [radian]')
plt.legend(['Euler', 'Euler-Cromer', 'ODE'])
plt.savefig("thetas.png")
plt.clf()

plt.plot(x[::3], omegas[0][::3], '.', x[1::3], omegas[1][1::3], '.', x[2::3], omegas[2][2::3], '.')
pos = plt.gca().get_position()
pos.x0 = 0.15
plt.gca().set_position(pos)
plt.title(f'Wartości prędkości kątowej dla metod')
plt.xlabel('t[s]')
plt.ylabel('Omega [radian/s]')
plt.legend(['Euler', 'Euler-Cromer', 'ODE'])
plt.savefig("omegas.png")
plt.clf()

# Zad 2
methods = ['Euler', 'Euler-Cromer', 'ODE']

# theta_diffs = []
thetas_np = np.array(thetas)
omegas_np = np.array(omegas)
for i in range(3):
    for j in range(i):
        # theta_diffs.append([[i,j],thetas_np[i]-thetas_np[j]])
        print(methods[i], methods[j], np.average((thetas_np[i]-thetas_np[j])**2), np.average((omegas_np[i]-omegas_np[j])**2))


# for indexes, diffs in theta_diffs:
#     plt.plot(list(range(steps)), diffs, label=f'{indexes[0]}-{indexes[1]}')
# plt.legend()
# plt.show()


# Zad 3

# [1E11, 1E10, 1E9, 1E8 ,1E7, 1E3, 1E2, 10, 5,4,3,2,1]

thetas = [pi/12*i for i in np.arange(0,1.01,0.01)]

# diffs = [[diff,{}] for diff in thetas]

EC_data = []
ODE_data = []

for theta in thetas:
    EC_sin = []
    EC_nosin = []
    ec = EulerCromer(0, theta)
    ec_nosin = EulerCromer(0, theta)
    for i in range(steps):
        EC_sin.append(ec.theta)
        EC_nosin.append(ec_nosin.theta)
        ec.iterate()
        ec_nosin.iterate_nosin()

    EC_data.append(np.average((np.array(EC_sin) - np.array(EC_nosin))**2))

    ode = ODE(0, theta)
    ODE_out = ode.get_sim(steps)
    ODE_t_sin = [ODE_out[i][1] for i in range(steps)]
    ode = ODE(0, theta)
    ODE_out = ode.get_sim_nosin(steps)
    ODE_t_nosin = [ODE_out[i][1] for i in range(steps)]
    
    ODE_data.append(np.average((np.array(ODE_t_sin) - np.array(ODE_t_nosin))**2))

plt.plot(thetas, EC_data, '.')
plt.plot(thetas, ODE_data, '.')
plt.title(f'Błąd między metodą z sinusem a bez')
plt.xlabel('Wstępne odchylenie [radian]')
plt.ylabel('Błąd')
plt.legend(['Euler-Cromer', 'ODE'])
plt.savefig('errors.png')

