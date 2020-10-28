import numpy as np
import matplotlib.pyplot as plt
import sdf
import os
import scipy.signal as signal
import constant as const
import function as func
import scipy.fftpack as fftpack
from matplotlib.ticker import MultipleLocator, FuncFormatter
from multiprocessing.dummy import Pool as ThreadPool
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import pyplot as plt

#import wigner
plt.switch_backend('agg')
###
####
def reg_cmap_transparent(iname,alpha):
	oname = iname + '_transparent'
	cmap = plt.get_cmap(iname)
	values = np.linspace(0,1,256)
	colors = cmap(values)
	for i in range(256):
		colors[i][3] = alpha[i]
	colorlist = [(values[i],colors[i]) for i in range(256)]
	cmap = plt.cm.colors.LinearSegmentedColormap.from_list(oname,colorlist)
	plt.cm.register_cmap(cmap=cmap)
	return cmap

def create_alpha(func):
	return [ 1 if func(i)>1 else 0 if func(i)<0 else func(i) for i in range(256)]

font = {'family' : 'monospace',  
	'color'  : 'black',  
	'weight' : 'normal',  
	'size'   : 15,  
	}  
color='hsv'
	index = 1
	plt.set_cmap(reg_cmap_transparent(color,create_alpha(lambda x:(1-abs(x/127.5-1)**index))))
cmap = plt.get_cmap()

	name='Ek_energe  '
	font_size = 20
fig = plt.figure(figsize=(12,8))
	ax = fig.add_subplot(111,projection='3d')
	ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([12,8,8,12]))
ax.view_init(45,-100)


def draw(x):
	savefigdir=const.figdir+str(x)
	sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
	data=sdf.read(sdfdir,dict=True)
	Bz=data['Electric Field/Ey']
	bz=Bz.data	
	density=data['Derived/Number_Density/electron1'].data
	bz=bz.T
	density=density
	#np.savetxt(const.txtdir+str(x)+'.txt',bz)	
	x    = data['Grid/Grid_mid'].data[0]/1.e-6
	y    = data['Grid/Grid_mid'].data[1]/1.e-6
	z    = data['Grid/Grid_mid'].data[2]/1.e-6
	#var1  = data['Magnetic Field/By_averaged'].data/bxunit
	#var2  = data['Magnetic Field/Bz_averaged'].data/bxunit
	#var   = (var1**2+var2**2)**0.5
	Y, X, Z= np.meshgrid(y, x, z)#, sparse=False, indexing='ij')
	print(X.shape)
	X.reshape(-1)
	Y.reshape(-1)
	Z.reshape(-1)
	density.reshape(-1)
	im=ax.scatter(X.reshape(-1),Y.reshape(-1),Z.reshape(-1),s=0.1,c=density.reshape(-1),cmap=cm.rainbow)
	cbar=plt.colorbar(im,ticks=np.linspace(450, 650, 5),fraction=0.015,pad=-0.1)
	cbar.set_label(name+r'$[meV]$',fontdict=font)
	fig.savefig('3db.png',dpi=400)

draw(2000)
