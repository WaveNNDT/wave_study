import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# %%
t = np.linspace(0, 100, 1000)  # time
k = 1
m = 1
omega_n = np.power(k / m, 0.5)  # natural angular frequency
c1 = 1
c2 = 1
c = 1   # damper
c_cr = 2 * omega_n * m

zeta = c / c_cr

# chi

u = c1 * np.exp((-zeta * omega_n + 1j * omega_n * np.power(1 - np.power(zeta, 2), 0.5) * t)) + c2 * np.exp(
    (-zeta * omega_n - 1j * omega_n * np.power(1 - np.power(zeta, 2), 0.5)) * t)

# ims = []
# for i in range(len(t)):

plt.figure()
plt.title(str(zeta))
plt.scatter(t, u)
plt.grid()

