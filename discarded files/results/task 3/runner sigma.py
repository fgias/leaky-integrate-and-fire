import numpy as np
import time
from matplotlib import pyplot as plt
from numpy import asarray
from numpy import savetxt
import csv
from notify_run import Notify

notify = Notify()


def run(N, R, sigma, time_max):
	mu = 1
	u_th = 0.98
	dt = 0.005
	
	seed_num = 567632
	np.random.seed(seed_num)
	u = np.random.rand(N)*u_th
	circles = np.zeros(N)

	t = 0
	num = 0
	plotfreq = 5000
	start_time = time.time()
	print('--------------------------')
	print('Starting calculation')
	print('--------------------------')
	while t < time_max:
		u = u + dt*(mu) # Modified model, lambda = 0.
		for i in range(N):
			coupling = 0 # R <= N/2.
			if i-R < 0:
				for j in range((i-R)%N,N):
					coupling = coupling + sigma*(u[i]-u[j])/(2*R)

				for j in range(i+R+1):
					coupling = coupling + sigma*(u[i]-u[j])/(2*R)

				u[i] = u[i] + dt*coupling
			elif i+R > N-1:
				for j in range(i-R,N):
					coupling = coupling + sigma*(u[i]-u[j])/(2*R)

				for j in range((i+R)%N+1):
					coupling = coupling + sigma*(u[i]-u[j])/(2*R)

				u[i] = u[i] + dt*coupling
			else: 
				for j in range(i-R,i+R+1):
					coupling = coupling + sigma*(u[i]-u[j])/(2*R)

				u[i] = u[i] + dt*coupling

			if u[i] > u_th:
				u[i] = 0
				circles[i] += 1

		if num%plotfreq == 0:
			print('t =', round(t))
			plt.scatter(range(1,N+1), u)
			plt.ylim([0, 1])
			plt.yticks([0, 0.5, 1])
			plt.ylabel("ui")
			plt.xticks([1, N/2, N])
			plt.xlabel("index i")
			plt.show()

		t = t + dt
		num = num + 1

	print('--------------------------')
	print('End of calculation')
	print('--------------------------')
	time_calc = time.time() - start_time
	print("Time: %s seconds." % (round(time_calc, 2)))
	print('--------------------------')

	"""
	Saving data.
	
	"""
	w = 2*np.pi*circles/t
	plt.scatter(range(1, N+1), w)
	plt.ylabel("ωi")
	plt.yticks([round(min(w),2), round((min(w)+max(w))/2,2), round(max(w),2)])
	plt.xticks([1, N/2, N])
	plt.xlabel("index i")
	plt.savefig('data/w_sigma=%s.png' % sigma, dpi=800)
	plt.show()

	plt.scatter(range(1,N+1), u)
	plt.ylim([0, 1])
	plt.yticks([0, 0.5, 1])
	plt.ylabel("ui")
	plt.xticks([1, N/2, N])
	plt.xlabel("index i")
	plt.savefig('data/u_sigma=%s.png' % sigma, dpi=800)
	plt.show()

	delta_w = max(w) - min(w)

	with open('data/deltaw_data.csv', 'a', newline='') as file:
		writer = csv.writer(file)
		writer.writerow([delta_w])

	u_data = asarray(u)
	w_data = asarray(w)
	savetxt('data/u_data_sigma=%s.csv' % sigma, u_data, delimiter=',')
	savetxt('data/w_data_sigma=%s.csv' % sigma, w_data, delimiter=',')
	with open('data/calc_sigma=%s.txt' % sigma, 'w') as f:
		f.write('Results: \n')
		f.write('delta_w = %s \n' % delta_w)
		f.write('time_calc = %s \n \n' % time_calc)
		f.write('Parameters: \n')
		f.write('N = %s \n' % N)
		f.write('R = %s \n' % R)
		f.write('sigma = %s \n' % sigma)
		f.write('time_max = %s \n' % time_max)
		f.write('seed_num = %s \n' % seed_num)

	notify.send('Calculation with sigma=%s ended.' % sigma)


for sigma in [1.9]:
	run(60, round(0.35*60), round(sigma,1), 3000)

