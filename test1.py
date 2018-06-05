import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# %%
t = np.linspace(0, 100, 1000)  # time
k = 0.5
m = 2
omega_n = np.power(k / m, 0.5)  # natural angular frequency
c1 = 2
c2 = 2

#%%
c = 0.1   # damper
c_cr = 2 * omega_n * m
zeta = c / c_cr

if np.isnan(np.power(1 - np.power(zeta, 2), 0.5)):
    u = c1 * np.exp((-zeta * omega_n + omega_n * np.power(np.power(zeta, 2) - 1, 0.5)) * t) + c2 * np.exp(
        (-zeta * omega_n + omega_n * np.power(np.power(zeta, 2) - 1, 0.5)) * t)
else:
    u = c1 * np.exp((-zeta * omega_n + 1j * omega_n * np.power(1 - np.power(zeta, 2), 0.5)) * t) + c2 * np.exp(
        (-zeta * omega_n - 1j * omega_n * np.power(1 - np.power(zeta, 2), 0.5)) * t)



# ims = []
# for i in range(len(t)):

f1 = plt.figure(1)
f1.clf()
f1.set_size_inches(5,5)
# plt.subplot(121)
plt.title(str(zeta))
plt.scatter(t, u, s=10)
plt.grid()
plt.xlabel('Time')
plt.ylabel('Amplitude')

# plt.subplot(122)
# plt.scatter(t, np.imag(u), s=10)
# plt.grid()


