# -- coding: utf-8 --
import os
#constant
nperseg = 256
name    = 'a2_n1'
data_name = "a2_n1/"
filenumber = 4
sdfdir  =  "../Data/"+data_name
txtdir  =  "txt/"+data_name
figdir  =  "fig/"+data_name
def checkdir():
	if (os.path.isdir(txtdir) == False):
    		os.mkdir(txtdir)
	if (os.path.isdir(figdir) == False):
    		os.mkdir(figdir)
checkdir()
###
c       =  3e8
micron  =  1e-6
lamada  =  0.8 * micron 
gridnumber = 900
Ny      =  250
Nx      =  gridnumber
Nz      =  200

start   =  1
stop    =  3
step    =  1
dt_snapshot= 2e-15
dt      =  dt_snapshot*1e15      #fs
x_max   =  36 * micron #lamada#* micron   #60 * lamada #micron
x_min   =  0  * lamada#* micron
x_end   =  x_max - x_min 
y_lenth =  20 * micron#lamada
z_lenth =  20 * micron
window_start_time =  (x_max - x_min) / c
delta_x =  x_end/gridnumber
t_end   =  stop * dt_snapshot
x_interval=1
t_total=1e15*x_end/c         #fs
t_size=t_total/(dt_snapshot*1e15)+1+1           #t_grid_number
######t_size=int(1e15*gridnumber*delta_x/c)+1

if t_end-window_start_time<0:
      xgrid   =  int(gridnumber)
else:
      xgrid   =  int(gridnumber + c*(t_end-window_start_time)/delta_x)
#####fft freqs
