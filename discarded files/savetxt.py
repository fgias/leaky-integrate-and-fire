from numpy import asarray
from numpy import savetxt
from numpy import loadtxt

udata = asarray(u)
timearray = np.array([t])
timedata = asarray(timearray)
circdata = asarray(circles)
savetxt('data/udata.csv', udata, delimiter=',')
savetxt('data/timedata.csv', timedata, delimiter=',')
savetxt('data/circledata.csv', circdata, delimiter=',')