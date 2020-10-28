# -- coding: utf-8 --
import math
import sdf
import math
import numpy as np
import os
import multiprocessing 
import time
from numpy import ma
import constant as const
import os
#from multiprocessing.dummy import Pool as ThreadPool 
savedir=const.txtdir   #"./txt/a0_1_2e-2/"
savename="xt.txt"
fftdir =const.figdir      #"./fig/a0_1_2e-2/"  
###
dirsdf  = const.sdfdir   # '../Data/a0_1_2e-2/'
dirsize =  const.filenumber    #4

bakdir="baktxt/"+const.data_name

b = const.x_max/3e8/const.dt_snapshot/2
b = int(b)
limit_min=0.1e12
limit_max=10e12

if (os.path.isdir(bakdir) == False):
        os.mkdir(bakdir)

def rm(n):
	save = 0
	if n==b:
		save = 1
	if n%10 == 0:
		save = 1
	if save == 0:
		if os.path.exists(dirsdf+str(n).zfill(dirsize)+".sdf")==True:
			print('begin_remove')
			#os.remove(dirsdf+str(n).zfill(dirsize)+".sdf")
			print('remove:')
	'''
	if save == 1:
		print('save_sdf:',i)
	'''	
	return
	#p "draw",x

def energe(x):
	limit_min=0.1e12
	limit_max=10e12
 
	#p "draw",x
	savefigdir=const.figdir+str(x)+'k_bz.png'
	sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
	data=sdf.read(sdfdir,dict=True)
	Bz=data['Electric Field/Ey']
	time=data['Header']['time']
	bz=Bz.data
	#bz=bz.T
	k_bz=np.fft.fft(bz,axis=0)
	delta_k=3.14/const.delta_x/(const.Nx/2)
	k_bz2=k_bz*1
	k_n=[]
	for n in range(0,const.Nx):
		mi = 3e8/limit_min
		ma = 3e8/limit_max
		if 2 * 3.14 / ma  > n * delta_k and  n * delta_k > 2 * 3.14 / mi:
			k_n.append(n)
	k_bz2[0:k_n[0],:,:]=0    #k_bz.argmin()
	k_bz2[k_n[-1]:-k_n[-1],:,:]=0  #k_bz.argmin()
	k_bz2[-k_n[0]:,:,:]=0    #k_bz.argmin()
	bz_filter=np.fft.ifft(k_bz2,axis=0)
	E_x=np.sum(np.sum(np.square(bz)))
	E_Thz=np.sum(np.sum(np.square(bz_filter.real)))
	E_2=np.array([E_x,E_Thz])
	np.savetxt("baktxt/"+const.data_name+str(x)+"E_2.txt",E_2)
	return
def baktxt(n):
	#print("save file:"+str(n))
	sdfdir=const.sdfdir +str(n).zfill(const.filenumber)+".sdf"
	data = sdf.read(sdfdir,dict=True)
	header=data['Header']
	time=header['time']
	Ex=data["Electric Field/Ex"].data
	Ex_y0=Ex[:,int(const.Ny/2),int(const.Nz/2)]
	Ey=data["Electric Field/Ey"].data
	Ey_y0=Ey[:,int(const.Ny/2),int(const.Nz/2)]
	ne=data['Derived/Number_Density/electron1'].data
	ne_y0=ne[:,int(const.Ny/2),int(const.Nz/2)]
	np.savetxt("baktxt/"+const.data_name+str(n)+"Ey_y0.txt",Ey_y0)
	np.savetxt("baktxt/"+const.data_name+str(n)+"Ex_y0.txt",Ex_y0)
	np.savetxt("baktxt/"+const.data_name+str(n)+"ne_y0.txt",ne_y0)
	return

def aaa(n):
	print('n',n)
	while 1:
		if os.path.exists(dirsdf+str(n+1).zfill(dirsize)+".sdf") == True:# or os.path.exists(dirsdf+str(stop).zfill(dirsize)+".sdf") == True:
	#time.sleep()
			break
	limit_min=0.1e12
	limit_max=10e12
	b = const.x_max/3e8/const.dt_snapshot/2
	b = int(b)
	sdfdir=const.sdfdir +str(n).zfill(const.filenumber)+".sdf"
	data=sdf.read(sdfdir,dict=True)
	energe(n)
	baktxt(n)
	rm(n)

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
	start   =  const.start
	stop    =  const.stop       #5889 #17000
	step    =  1
	dt_snapshot= const.dt_snapshot     #9e-15
	dt      =  dt_snapshot*1e15      #fs
	x_max   =  const.x_max      #80 * lamada   #60 * lamada    #micron
	x_min   =  0 * micron
	x_end   =  x_max - x_min
	y       =  const.Ny/2 
	z       =  const.Nz/2
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

	pool = multiprocessing.Pool(processes=32)
	results = pool.map(aaa,range(int(start),int(stop+step),int(step)))
	pool.close()
	pool.join()
	print('finish')


	#time = range(0,stop+step,step)
	#locate = (time*const.dt_snapshot - const.window_start_time)*3e8*1e6
	#np.savetxt(const.txtdir + 'eff_locate.txt',locate)
