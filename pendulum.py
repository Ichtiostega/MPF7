from math import pi, sin
from scipy.integrate import odeint
import numpy as np
from tabulate import tabulate
from matplotlib import pyplot as plt

dt = 0.1
g = 10
l = 100

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

steps = 2000

ode = ODE(0, pi/10000)
ode_out = ode.get_sim(steps)

data = []
e = Euler(0, pi/10000)
ec = EulerCromer(0, pi/10000)
print(e, ec)
for i in range(steps):
    data.append([e.omega, e.theta, ec.omega, ec.theta, ode_out[i][0], ode_out[i][1]])
    e.iterate()
    ec.iterate()

print(tabulate(data, headers=['e_o', 'e_t', 'ec_o', 'ec_t', 'ode_o', 'ode_t']))

x = [i*dt for i in range(0, steps)]

thetas = list(zip(*data))[1::2]
omegas = list(zip(*data))[::2]

plt.plot(x, thetas[0], x, thetas[1], x, thetas[2])
plt.show()

plt.plot(x, omegas[0], x, omegas[1], x, omegas[2])
plt.show()