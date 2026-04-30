import matplotlib.pyplot as plt
import numpy as np
import sxs

# Remember that it is always best to specify the version of the catalog you want to load:
df = sxs.load("dataframe", tag="3.0.0")

sxs_bhns = sxs.load("SXS:BHNS:0009")
w = sxs_bhns.h


w_2_2 = w[:, w.index(2, 2)]
w_3_3 = w[:, w.index(3, 3)]
t22 = w_2_2.t[np.argmax(w_2_2.abs)]

ratio = w_3_3.abs / w_2_2.abs

'''
plt.plot(w_2_2.t, w_2_2.abs)
plt.plot(w_2_2.t, w_2_2)
plt.plot(w_3_3.t, w_3_3.abs)
plt.plot(w_3_3.t, w_3_3)
'''
plt.plot(w_2_2.t+10, ratio)
plt.xlim([t22,w_2_2.t[-1]])
plt.xlabel(r"$(t_{\mathrm{corr}} - r_\ast)/M$")
plt.ylabel(r"$h^{(3,3)} / h^{(2,2)}$")
plt.show()
