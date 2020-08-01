# -- coding: utf-8 --
import math
import sdf
import matplotlib
import math
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import multiprocessing 
import time
from numpy import ma
from matplotlib import colors, ticker, cm
#from matplotlib.mlab import bivariate_normal
#from scipy.interpolate import spline
import constant as const
from multiprocessing import shared_memory
#from multiprocessing import Process, Array,Value
#from multiprocessing.dummy import Pool as ThreadPool 
savedir=const.txtdir   #"./txt/a0_1_2e-2/"
savename="xt.txt"
fftdir =const.figdir      #"./fig/a0_1_2e-2/"  
###
dirsdf  = const.sdfdir   # '../Data/a0_1_2e-2/'
dirsize =  const.filenumber    #4
def extract(n):
        #### header data ####
	print('n:'+str(n))
	#xt = np.frombuffer(global_arr_shared, np.double).reshape(SHAPE)
	#global a.shape
	#global a.dtype
	xt = np.ndarray(ss1, dtype=ss2, buffer=shm.buf)
	data = sdf.read(dirsdf+str(n).zfill(dirsize)+".sdf",dict=True)
	header=data['Header']
	time=header['time']
	E_y0=data['Electric Field/Ey'].data[:,int(y)]
	if  n  <  start_move_number:
		for x in range(1,int(gridnumber/x_interval)+1):
			a=int(x*x_interval)
			d_n=int((1e15*delta_x*a/c)/dt)
			if n-d_n > 0 and n-d_n < t_size :
#[fs]
				xt[x][n-d_n]=E_y0[a-1] #/bxunit            
	else:
		for x in range(1,int(xgrid/x_interval)+1):

			a=int(x*x_interval)
			if a-c*(time-window_start_time)/delta_x >= 0 and a-c*(time-window_start_time)/delta_x < gridnumber-1:
	#[fs]
				d_n=int((1e15*delta_x*a/c)/dt)
				xt[x][n-d_n]=E_y0[int(round(a-c*(time-window_start_time)/delta_x))]  #/bxunit
	return "OK"+str(n)
                   #else:bz.append(0)
                   #print 'Reading finished%d' %len(t)
if __name__ == "__main__":
	######## Constant defined here ########
	pi        =     3.1415926535897932384626
	q0        =     1.602176565e-19 # C
	m0        =     9.10938291e-31  # kg
	v0        =     2.99792458e8    # m/s^2
	kb        =     1.3806488e-23   # J/K
	mu0       =     4.0e-7*pi       # N/A^2
	epsilon0  =     8.8541878176203899e-12 # F/m
	h_planck  =     6.62606957e-34  # J s
	####lamada


	wavelength=     const.lamada     #10.6e-6

	####

	frequency =     v0*2*pi/wavelength
	micron    =     1e-6
	c         =     3e8
	exunit    =     m0*v0*frequency/q0
	bxunit    =     m0*frequency/q0
	denunit    =     frequency**2*epsilon0*m0/q0**2
	print('electric field unit: '+str(exunit))
	print('magnetic field unit: '+str(bxunit))
	print('density unit nc: '+str(denunit))
	font = {'family' : 'monospace',  
		'color'  : 'black',  
		'weight' : 'normal',  
		'size'   : 28,  
	}  
	if (os.path.isdir(savedir) == False):
		os.mkdir(savedir)
		
	if (os.path.isdir(fftdir) == False):
		os.mkdir(fftdir)    
	######### Script code drawing figure ################
	######constant
	###
	c       =  3e8
	micron  =  1e-6
	lamada  =  const.lamada #10.6 * micron
	gridnumber = const.Nx     #2400
	start   =  1
	stop    =  const.stop       #5889 #17000
	step    =  1
	dt_snapshot= const.dt_snapshot     #9e-15
	dt      =  dt_snapshot*1e15      #fs
	x_max   =  const.x_max      #80 * lamada   #60 * lamada    #micron
	x_min   =  0 * micron
	x_end   =  x_max - x_min
	y       =  const.Ny/2 
	window_start_time =  (x_max - x_min) / c
	#start_move_number = window_start_time * 1e15      #fs
	start_move_number =  int(window_start_time / dt_snapshot)
	delta_x =  x_end/gridnumber
	t_end   =  stop * dt_snapshot
	x_interval=const.x_interval          #10
	t_total=1e15*x_end/c         #fs
	t_size=t_total/(dt_snapshot*1e15)+1           #t_grid_number
	if t_end-window_start_time<0:
		xgrid   =  int(gridnumber)
	else: 
		xgrid   =  int(gridnumber + c*(t_end-window_start_time)/delta_x)

####################
	x_interval=const.x_interval      #10
	t_total=1e15*x_end/c         #fs
	t_size=int(t_total/dt)+1+1   

######allay define
	SHAPE = ((int(xgrid/x_interval)+1,t_size))
	#xt = Array('f',SHAPE)
	#global a
	a=np.zeros((int(xgrid/x_interval)+1,t_size))
	shm = shared_memory.SharedMemory(create=True, size=a.nbytes)
	#xt = np.ndarray(a.shape, dtype=a.dtype, buffer=shm.buf)
	#xt[:,:]=a[:,:]
	#import multiprocessing 
	#import time
	#import numpy as np 


	#global_arr_shared = None

	#def init_pool(arr_shared):
	#	global global_arr_shared
	#	global_arr_shared = arr_shared

	#def worker(i):
	#	arr = np.frombuffer(global_arr_shared, np.double).reshape(SHAPE)
	#	time.sleep(1)  # some other operations
	#	return np.sum(arr * i)

	#arr = np.zeros(SHAPE)
	#arr_shared = multiprocessing.RawArray('d', arr.ravel())
	#with multiprocessing.Pool(processes=48, initializer=init_pool, initargs=(arr_shared,)) as pool:  # initargs传入tuple
	#	for result in pool.map(extract,range(start,stop+step,step)):
	#		print(result)
  #pool = ThreadPool(48)
	ss1,ss2=a.shape,a.dtype
	pool = multiprocessing.Pool(processes=4,initargs=(ss1,ss2))
	#for i in range(start,stop+step,step):
	#	results.append(pool.apply_async(extract, (i, ))) 
	results = pool.map(extract,range(int(start),int(stop+step),int(step)))

	pool.close()
	pool.join()

  #for n in range(start,stop+step,step):
	#extract(n)
		#for res in results:
			#print(res.get())
	#arr = np.array(xt)
	xt = np.ndarray(a.shape, dtype=a.dtype, buffer=shm.buf)
	np.savetxt(savedir+savename, xt)
	shm.close()
	shm.unlink()
