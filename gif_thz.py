import numpy as np
import matplotlib.pyplot as plt
import sdf
import os
import scipy.signal as signal
import constant as const
import imageio
import function as func
import scipy.fftpack as fftpack
from matplotlib.ticker import MultipleLocator, FuncFormatter
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing
#import wigner
plt.switch_backend('agg')
###
####
interval=10
image_list=[]
png_savedir = "./gif/png/"

def x_formatter(x, pos):
        a=(const.delta_x*x + const.c*T)   *1e6 
        return  "%d"%int(a)
def x_formatter2(x, pos):
	#locate = x * const.nperseg
        a=x + const.c*T * 1e6
        return  "%d"%int(a)

def draw(x):
	#savefigdir=const.figdir+str(int(locate/1e-6))+'_'+str(Thz1)+'_'+str(Thz2)+'k_bz.png'
	sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
	data=sdf.read(sdfdir,dict=True)
	Bz=data['Electric Field/Ey']
	time=data['Header']['time']
###
	global T
	time=data['Header']['time']
	if time-const.window_start_time<0:
		T=0
	else:
		T=time-const.window_start_time
###



	bz=Bz.data[:,:,int(const.Nz/2)]
	k_bz=np.fft.fft(bz,axis=0)
	k_bz2d=np.fft.fft2(bz)
	delta_k=3.14/const.delta_x/(const.Nx/2)
	k_bz2=k_bz*1
	k_n=[]
	for n in range(0,const.Nx):
		mi = 3e8/0.1e12
		ma = 3e8/10e12
		if 2 * 3.14 / ma  > n * delta_k and  n * delta_k > 2 * 3.14 / mi:
			k_n.append(n)
	print("n",k_n[0],k_n[-1])
	k_bz2[0:k_n[0],...]=0    #k_bz.argmin()
	k_bz2[k_n[-1]:-k_n[-1],...]=0  #k_bz.argmin()
	k_bz2[-k_n[0]:,...]=0    #k_bz.argmin()

	fs=2*3.14/const.delta_x

	bz_filter=np.fft.ifft(k_bz2,axis=0)
	E_x=np.sum(np.sum(np.square(bz)))
	E_Thz=np.sum(np.sum(np.square(bz_filter.real)))
	eff=E_Thz/E_x
	print("efficiency",E_x,E_Thz,eff)
	ne=data['Derived/Number_Density/electron1'].data[:,:,int(const.Nz/2)]
	ne_y0=ne[:,int(const.Ny/2)]
	ne=ne.T	
	fig,axs=plt.subplots(2,2)
	im=axs[0][0].pcolormesh(bz.T,cmap=plt.get_cmap('bwr'))   
	im2=axs[1][0].pcolormesh(ne)
	im3=axs[0][1].pcolormesh(bz_filter.real.T,cmap=plt.get_cmap('bwr'))
	fig.colorbar(im,ax=axs[0][0])
	fig.colorbar(im3,ax=axs[0][1])
	im4=axs[1][1].pcolormesh(abs(k_bz2d[:int(const.Nx/2),:int(const.Ny/2)].T),cmap=plt.get_cmap('bwr'))
	#fig.savefig(savefigdir,dpi=200)

	t_fs=int(time/1e-15)
	axs[0][0].set_title("t="+str(t_fs),fontsize=12,color='r')
	#axs[0][0].set_title(str(x*const.dt_snapshot))
	axs[0][0].xaxis.set_major_formatter( FuncFormatter( x_formatter ) )

	fig.savefig(png_savedir+str(x)+"ref_k.png",dpi=200)
	image_list.append(png_savedir+str(x)+"ref_k.png")
	plt.clf()
	plt.close('all')


#wigner.aaa() 
'''        
for i in range(0,const.stop,interval):
        draw(i)
'''
image_list=[]
'''
pool = multiprocessing.Pool(processes=4)   #,initargs=(ss1,ss2))
               #for i in range(start,stop+step,step):
        #       results.append(pool.apply_async(extract, (i, ))) 
results = pool.map(draw,range(0,const.stop,interval))
'''
results=[]
for i in range(0,const.stop,interval):
	results.append(draw(i))

for x in range(0,const.stop,interval):
        image_list.append(png_savedir+str(x)+"ref_k.png")
'''
piool.close()
pool.join()
'''


def create_gif(image_list, gif_name, duration = 0.5):
    '''
    '''
    frames = []
    #for image_name in image_list:
        #frames.append(imageio.imread(image_name))

    #imageio.mimsave(gif_name, frames, 'GIF', duration=duration)
 
    pool = ThreadPool()
    frames = pool.map(imageio.imread,image_list)
    imageio.mimsave(gif_name, frames, 'GIF', duration=duration)
    pool.close()
    pool.join()
    return

def main():
    #image_list = ['1.jpg', '2.jpg', '3.jpg']
    gif_name = "gif/"+const.name+"_thz.gif"
    duration = 0.1
    create_gif(image_list, gif_name, duration)

#if __name__ == '__main__':
const.checkdir()
main()
