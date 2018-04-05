import matplotlib.pyplot as plt
import numpy as np
import math

#Loading data 
rhoInitial = np.loadtxt("rhoInitial.txt")
rhoFinal = np.loadtxt("rhoFinal.txt")
pressure=np.loadtxt("pressureFinal.txt")

# Figure 1
plt.figure(1)
plt.imshow(rhoInitial)
plt.colorbar()
plt.savefig("rhoInitial.png")
plt.show()

# Figure 2
plt.figure(2)
plt.imshow(rhoFinal)
plt.savefig("rhoFinal.png")
plt.colorbar()
plt.show()

