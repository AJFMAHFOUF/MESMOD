import numpy as np
import matplotlib.pyplot as plt
i,j,x=np.loadtxt('../data/fort.112',unpack=True)
y=np.reshape(x,(26,26),order='F')
plt.contourf(y)
plt.colorbar()
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('Initial field')
plt.show()
