#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

N   = 128
pi  = np.pi
x   = np.linspace(0, 2.0*pi, N)
y   = np.sin(pi*x)

plt.plot(x,y,'r', lw=1.5)
plt.xlabel('$0 \leq x \leq 2$')
plt.ylabel('$f(x)$')
print("Use this script to test keys")
plt.show()
