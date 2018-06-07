import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# %%
t = np.linspace(0, 50, 1000)  # time
k = 1
m = 2
omega_n = np.power(k / m, 0.5)  # natural angular frequency

c1 = 0.5
c2 = 0.5

c = [0.5, 1, 2, 5]  # damper
c_cr = 2 * omega_n * m
zeta = c / c_cr

u = np.empty([len(c), len(t)], np.complex128)
for i in range(len(c)):
    omega_d = omega_n * np.power(1 - np.power(zeta, 2), 0.5)  # damped natural frequency
    if np.isnan(omega_d[i]):
        omega_d[i] = (omega_n * np.power(np.power(zeta[i], 2) - 1, 0.5))
        u[i] = c1 * np.exp((-zeta[i] * omega_n + omega_d[i]) * t) + c2 * np.exp((-zeta[i] * omega_n - omega_d[i]) * t)
        # b = np.exp(-zeta * t) * np.cos(omega_d * t)
    else:
        u[i] = c1 * np.exp((-zeta[i] * omega_n + 1j * omega_d[i]) * t) + c2 * np.exp(
            (-zeta[i] * omega_n - 1j * omega_d[i]) * t)
        # b = np.exp(-zeta * t) * np.cos(omega_d * t)

f1 = plt.figure(1)
f1.clf()
f1.set_size_inches(10, 5)
plt.plot(t, u[0], label='zeta = ' + str(np.round(zeta[0], 2)))
plt.plot(t, u[1], label='zeta = ' + str(np.round(zeta[1], 2)))
plt.plot(t, u[2], label='zeta = ' + str(np.round(zeta[2], 2)))
plt.plot(t, u[3], label='zeta = ' + str(np.round(zeta[3], 2)))
# plt.scatter(t, u[0], s=10)
# plt.scatter(t, u[1], s=10)
# plt.scatter(t, u[2], s=10)
# plt.scatter(t, u[3], s=10)
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.grid()
plt.legend(loc=0)
f1.tight_layout()

plt.savefig('Figures/comparison.png', dpi=600)
#
# F_hat = 1
# omega = 2
# F = F_hat * np.cos(omega * t)
# f_hat = F_hat / m
#
# # %%
# if np.isnan(np.power(1 - np.power(zeta, 2), 0.5)):
#     u = c1 * np.exp((-zeta * omega_n + omega_n * np.power(np.power(zeta, 2) - 1, 0.5)) * t) + c2 * np.exp(
#         (-zeta * omega_n + omega_n * np.power(np.power(zeta, 2) - 1, 0.5)) * t)
# else:
#     u = c1 * np.exp((-zeta * omega_n + 1j * omega_n * np.power(1 - np.power(zeta, 2), 0.5)) * t) + c2 * np.exp(
#         (-zeta * omega_n - 1j * omega_n * np.power(1 - np.power(zeta, 2), 0.5)) * t)
#
# # ims = []
# # for i in range(len(t)):
#
#
# # p = (1/(-np.power(omega,2)+1j*2*zeta*omega_n*omega+np.power(omega_n,2)))*f_hat*np.exp(1j*omega*t)
# p = (1 / (-np.power(omega, 2) + 1j * 2 * zeta * omega_n * omega + np.power(omega_n, 2))) * f_hat * np.cos(omega * t)
# p2 = (1 / (-np.power(omega, 2) + 1j * 2 * zeta * omega_n * omega + np.power(omega_n, 2))) * f_hat * np.exp(
#     1j * omega * t)
#
# f1 = plt.figure(1)
# f1.clf()
# f1.set_size_inches(10, 5)
# plt.subplot(121)
# plt.title(str(zeta))
# # plt.scatter(t, u, s=10)
# plt.scatter(t, u, s=10)
# plt.xlabel('Time')
# plt.ylabel('Amplitude')
# plt.grid()
#
# plt.subplot(122)
# plt.scatter(t, v, s=10)
# plt.grid()
# plt.xlabel('Time')
# plt.ylabel('Amplitude')
#
# # plt.subplot(122)
# # plt.scatter(t, np.imag(u), s=10)
# # plt.grid()
#
# # %%
# f2 = plt.figure(2)
# f2.clf()
# plt.scatter(t, b + np.absolute(p), s=10)
# plt.scatter(t, b + np.absolute(p2), s=10)
#
# p3 = p2 - p
# plt.figure()
# plt.scatter(t, np.real(p3))
# plt.scatter(t, np.imag(p3))
