import numpy as np
import matplotlib.pyplot as plt

kpi = 0.0008
h0 = 1 * np.pi
p = np.linspace(0.12, 2*np.pi - 0.12, 1000)
pr = np.zeros_like(p)
pk = np.zeros_like(p)
pe = np.zeros_like(p)

def Spring(p, h0, kpi, _):
    r = 5
    k = 5
    e = 5
    return r, k, e

for i in range(len(p)):
    r, k, e = Spring(p[i], h0, kpi, 1)
    pr[i] = r
    pk[i] = k
    pe[i] = e

plt.figure()
plt.plot(p, pr)
ax = plt.gca()
# ax.set_fontsize(14)
ax.set_xticks([0, 0.5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
ax.set_xticklabels(['0', '0.5π', 'π', '1.5π', '2π'])
# plt.ylim([-np.inf, np.inf])
plt.xlim([0, 2*np.pi])
plt.grid(True)