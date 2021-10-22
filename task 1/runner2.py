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
	dt = 0.01
	
	seed_num = 38958
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
		for i in range(N):
			u[i] = u[i] + dt*(mu - u[i])
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
		f.write('Results:')
		f.write('\n')
		f.write('delta_w = %s' % delta_w)
		f.write('\n')
		f.write('time_calc = %s' % time_calc)
		f.write('\n')
		f.write('\n')
		f.write('Parameters:')
		f.write('\n')
		f.write('N = %s' % N)
		f.write('\n')
		f.write('R = %s' % R)
		f.write('\n')
		f.write('sigma = %s' % sigma)
		f.write('\n')
		f.write('time_max = %s' % time_max)
		f.write('\n')
		f.write('seed_num = %s' % seed_num)

	notify.send('Calculation with sigma=%s ended.' % sigma)

N=70
for sigma in np.arange(0.2, 3.1, 0.2):
	run(N, round(0.4*N), round(sigma,1), 1000)


