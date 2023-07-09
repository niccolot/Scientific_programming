from matplotlib import pyplot as plt
import numpy as np
from itertools import groupby
from scipy.optimize import curve_fit

def fit_func(x,a):
    return x**a

def get_data(file):

    list = []
    with open(file) as f:
        for k, g in groupby(f, lambda x: x.startswith('\n')):
            if not k:
                list.append(np.array([[float(x) for x in d.split()] for d in g if len(d.strip())]))

    return np.mean(np.array(list), axis=0)


n_1000 = get_data('n1000.txt')
n_10000 = get_data('n10000.txt')
n_100000 = get_data('n100000.txt')
n_1000000 = get_data('n1000000.txt')

Smean_c1 = [np.squeeze(n_1000[np.where(n_1000[:,0]==1),1]), 
           np.squeeze(n_10000[np.where(n_10000[:,0]==1),1]),
           np.squeeze(n_100000[np.where(n_10000[:,0]==1),1]),
           np.squeeze(n_1000000[np.where(n_1000000[:,0]==1),1])]

Smax_c1 = [np.squeeze(n_1000[np.where(n_1000[:,0]==1),2]), 
           np.squeeze(n_10000[np.where(n_10000[:,0]==1),2]),
           np.squeeze(n_100000[np.where(n_10000[:,0]==1),2]),
           np.squeeze(n_1000000[np.where(n_1000000[:,0]==1),2])]

N_values = [1000,10000,100000,1000000]

pars_smean, _ = curve_fit(fit_func, N_values, Smean_c1)
pars_smax, _ = curve_fit(fit_func, N_values, Smax_c1)
x = np.linspace(1000,1000000,100000)


plt.figure(figsize=(10,10))

plt.subplot(2,2,1)
plt.plot(n_1000[:,0], n_1000[:,1], '.', label='N 1000')
plt.plot(n_10000[:,0], n_10000[:,1], '.', label='N 10000')
plt.plot(n_100000[:,0], n_100000[:,1], '.', label='N 100000')
plt.plot(n_1000000[:,0], n_1000000[:,1], '.', label='N 1000000')
plt.xlabel('c')
plt.ylabel("$<S'_{mean}>$")
plt.yscale('log')
plt.legend()

plt.subplot(2,2,2)
plt.plot(n_1000[:,0], n_1000[:,2], '.', label='N 1000')
plt.plot(n_10000[:,0], n_10000[:,2], '.', label='N 10000')
plt.plot(n_100000[:,0], n_100000[:,2], '.', label='N 100000')
plt.plot(n_1000000[:,0], n_1000000[:,2], '.', label='N 1000000')
plt.xlabel('c')
plt.ylabel("$<S_{max}>$")
plt.legend()

plt.subplot(2,2,3)
plt.plot(N_values,Smean_c1, '.')
plt.plot(x, fit_func(x, *pars_smean))
plt.xlabel('N')
plt.ylabel("$<S'_{mean}>_{c=1}$")
plt.xscale('log')
plt.yscale('log')

plt.subplot(2,2,4)
plt.plot(N_values,Smax_c1, '.')
plt.plot(x, fit_func(x, *pars_smax))
plt.xlabel('N')
plt.ylabel("$<S'_{max}>_{c=1}$")
plt.xscale('log')
plt.yscale('log')

plt.show()

