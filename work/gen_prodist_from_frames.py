#!/usr/local/bin/python
import numpy as np
import matplotlib.pyplot as plt

pot  = np.loadtxt("pot_300K.dat")
pot2 = np.loadtxt("pot_600K.dat")

plt.hist(pot, 50)
plt.hist(pot2, 50)
plt.xlim([-113100.000,-81800.000])
plt.ylim([1,500])

plt.show()
