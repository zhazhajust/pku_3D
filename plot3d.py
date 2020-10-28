import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import pyplot as plt
import sdf
plt.switch_backend('Agg')
#from matplotlib import animation
n=3200
from_path='../Data/a2_n1_T6_w8/'
to_path  ='./fig/'
def rebin3d(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1],shape[2],a.shape[2]//shape[2]
    return a.reshape(sh).mean(-1).mean(3).mean(1)


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
color='bwr'
index = 6
plt.set_cmap(reg_cmap_transparent(color,create_alpha(lambda x:(1-abs(x/127.5-1)**index))))
cmap = plt.get_cmap()

name='Ek_energe  '
font_size = 20
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111,projection='3d')
ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([12,8,8,12]))
ax.view_init(45,-90)
data=sdf.read(from_path+str(n).zfill(4)+'.sdf',dict=True)
header=data['Header']
time =header['time']
x    = data['Grid/Grid_mid'].data[0]/1.e-6
y    = data['Grid/Grid_mid'].data[1]/1.e-6
z    = data['Grid/Grid_mid'].data[2]/1.e-6
X, Y, Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')
color= data['Electric Field/Ey'].data
color  = color.reshape(np.size(color))
X    = X.reshape(np.size(X))
Y    = Y.reshape(np.size(Y))
Z    = Z.reshape(np.size(Z))

im=ax.scatter(X,Y,Z,s=0.1,c=color,cmap=cmap)

cbar=plt.colorbar(im,ticks=np.linspace(450, 650, 5),fraction=0.015,pad=-0.1)
cbar.set_label(name+r'$[meV]$',fontdict=font)
ax.set_xlabel('\n\nX'+ '[$\mu m$]',fontdict=font)
ax.set_ylabel('\n\ny'+ '[$\mu m$]',fontdict=font)
ax.set_zlabel('\n\nZ'+ '[$\mu m$]',fontdict=font)
fig.savefig('3dthz.png',dpi=400)
#def rotate(angle): 
#   ax.view_init(azim=angle) 
#rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0,362,2),interval=100) 
#rot_animation.save('rotation.gif', dpi=80, writer='imagemagick') 
