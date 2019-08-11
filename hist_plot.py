import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.stats import norm
import matplotlib.mlab as mlab

ifile = open(sys.argv[1],"r")

gap = []
for i in ifile:
    gap.append(i.strip())

gap = np.array(gap)

#x = np.asarray(np.linspace(float(min(gap)),float(max(gap)), num=int(len(gap))))                                                                                                                            

# Creating an array of data                                                                                                                                                                                 
x = np.asarray(list(map(lambda gap:gap.strip(),gap)), dtype='f')

#fitting best distribution to the data                                                                                                                                                                      
(mu, sigma) = norm.fit(x)

# Histogram                                                                                                                                                                                                 
n, bins, patches = plt.hist(x, bins='auto', normed=1, histtype='bar')

y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=2)

#plot the results                                                                                                                                                                                           
plt.xlabel('Homo-Lumo(eV)')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.grid(True)

plt.show()
