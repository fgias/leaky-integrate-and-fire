import numpy as np
import time
from matplotlib import pyplot as plt
from numpy import asarray
from numpy import savetxt
from numpy import loadtxt
import csv

t=2000
N=150

for lmd in np.arange(1,0.89,-0.01):
    w = loadtxt('w_data_lambda=%s_t=%s.csv' % (round(lmd,2), t), delimiter=',')
    plt.scatter(range(1, N+1), w)
    plt.ylabel("ωi")
    plt.yticks([2.65, 3.4])
    plt.xticks([1, N/2, N])
    plt.xlabel("index i")
    plt.savefig('for_report/w_lambda=%s_t=%s.png' % (round(lmd,2), round(t)), dpi=800)
    plt.close()
    
    delta_w = max(w) - min(w)
    with open('for_report/deltaw_data.csv', 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([delta_w])
    	