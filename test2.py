import numpy as np
import matplotlib.pyplot as plt

# %%
x = np.linspace(0, 1, 1000)
y = np.sin(x)
A1 = 1
A2 = 1.5
L = 1
t = 2
f1 = 1
f2 = 2

y1 = 0.5*A1*(np.cos((2*np.pi*x/L)-2*np.pi*f1*t) - np.cos((2*np.pi*x/L)+2*np.pi*f1*t))

plt.figure()
plt.plot(x, y1)

