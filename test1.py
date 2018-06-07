import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# %%
t = np.linspace(0, 100, 1000)  # time

#%% test
zeta = []
aa = [0.1, 1, 2, 3]
for i in range(len(aa)):
    zeta.append(aa[i] / 0.5)

u = np.empty([len(aa), len(t)], dtype=np.complex128)

for i in range(len(aa)):
    u[i] = aa[i] * np.exp(1*t+1)


#%%
plt.figure()
plt.clf()
plt.plot(t, u[0])
plt.plot(t, u[1])