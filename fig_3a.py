# -*- coding: UTF-8 -*-
####import sdf data#####
import constant as const
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

import sdf
import scipy.fftpack as fftpack
import numpy as np

n=100

sdfdir=const.sdfdir + str(n).zfill(4) + '.sdf'
data=sdf.read(sdfdir,dict=True)
bz=data['Electric Field/Ey'].data[:,:,int(188/2)]
ne=data['Derived/Number_Density/electron1'].data[:,:,int(188/2)]
ne_y0=ne[:,int(188/2)]
bz_y0=bz[:,int(188/2)]
hx = fftpack.hilbert(bz_y0)
hy = np.sqrt(bz_y0**2+hx**2)
import numpy
import sdf
import constant as const
import scipy.signal as signal
import numpy as np
pi=3.14



savedir = const.figdir + str(n)+'_fig2a.png'



def E_x_y_zxx(a):        
    #sdfdir=const.sdfdir +str(a).zfill(const.filenumber)+".sdf"
    #data=sdf.read(sdfdir,dict=True)
    Ex=data["Electric Field/Ex"].data[:,:,int(188/2)]
    Ex_y0=Ex[:,int(188/2)]
    Ey=data["Electric Field/Ey"].data[:,:,int(188/2)]
    Ey_y0=Ey[:,int(188/2)]
    k,x,zxx=signal.stft(Ey_y0,fs=2*pi/const.delta_x,nperseg=const.nperseg)
    zxx=abs(zxx)
    index = np.unravel_index(zxx.argmax(),zxx.shape)
    k_x=np.ones(const.Nx)*k[index[0]]
    for i in range(0,const.Nx):
        x=i
        a=zxx[...,int(x/const.nperseg)]
        a=a.tolist()
        a[0]=0
        max_index=a.index(max(a))
        if max(a) > 0.2 * zxx[index[0]][index[1]]:
            k_x[i]=(k[max_index])
    ne=data['Derived/Number_Density/electron1'].data
    ne_y0=ne[:,int(const.Ny/2),int(const.Ny/2)]

    return [Ex_y0,Ey_y0,zxx,ne_y0,k_x]
def k(x,zxx):
    a=zxx[...,int(x/const.nperseg)]
    a=a.tolist()
    a[0]=0
    max_index=a.index(max(a))
    # index=(max_index*pi)/(len(a)*const.delta_x)
    index =  (pi/const.delta_x)/len(a) * max_index
    #print "max(k_x),max_index,len(k_x),index",max(a),max_index,len(a),index # ,index
    return index
def scalar_p(Ex_y0):
    scalar_p=0
    e=1.6e-19
    c=3e8
    m0=9.1e-31
    import constant as const
    a=np.zeros(const.Nx)
    x=const.Nx - 1
    for i in range(0,x+1):
        scalar_p=Ex_y0[x-i]*const.delta_x+scalar_p 
    #a.append(e*scalar_p/(m0*c**2))   
        a[x-i]=e*scalar_p/(m0*c**2)  
    #   print  'n_scalar',a
    return a
def wp_2(x,ne_y0):
    ne=ne_y0[x]
    #ne=1e25
    e=1.6e-19
    m0=9.1e-31
#print "ne",ne
    w_p_2=ne*e**2/(m0*8.85e-12)
#print "w_p",w_p_2
    return w_p_2
def ref_index(x,Ex_y0,Ey_y0,zxx,ne_y0,k_x):
    c=3e8
    wp_x_2=wp_2(x,ne_y0)
#print "wp^2",wp_x_2
#p_y0=scalar_p(Ex_y0)
    p_y0_x=p_y0[x]
#print "scalar",p_y0_x
#k_x=k(x,zxx)
    k_x=k_x[x]
    w0=k_x*c/1
#print "w0",w0
    a=wp_x_2/(1+p_y0_x)
    b=(k_x*c)**2
    c=-(k_x*c)**2
    #print("a,b,c",a,b,c)
    result=np.roots([a,b,c])
    #print(result)
    return np.roots([a,b,c])
def ref_a(a):   
    ref=[]
    Ex_y0,Ey_y0,zxx,ne_y0,k_x=E_x_y_zxx(a)
    x=const.Nx
    global p_y0
    p_y0=scalar_p(Ex_y0)
    for i in range(0,1000):
        index=ref_index(i,Ex_y0,Ey_y0,zxx,ne_y0,k_x)
        #print("c",c)
        s_a= type(index) == numpy.ndarray
        #print("type",a)
        if type(index) == numpy.ndarray:
            index=index[-1]
        #print("c[-1]",c)
        ref.append(index)
    #print("ref:",type(ref[0]))
    return ref
#def help():
#print "func.ref_a(a)"
#print "a,x=1,1"
#print "Ex_y0,Ey_y0,zxx,ne_y0,k_x=func.E_x_y_zxx(a)"
#print "global p_y0"
#print "func.p_y0=scalar_p(Ex_y0)"    
#print "func.ref_index(x,Ex_y0,Ey_y0,zxx,ne_y0)"
c=3e8
'''
x=2
Ex_y0,Ey_y0,zxx,ne_y0,k_x=E_x_y_zxx(2800)
print(zxx)
wp_x_2=wp_2(x,ne_y0)
#print(wp_x_2)
#print "wp^2",wp_x_2
#p_y0=scalar_p(Ex_y0)
p_y0_x=p_y0[x]
#print(p_y0_x)
#print "scalar",p_y0_x
#k_x=k(x,zxx)
k_x=k_x[x]
#print(k_x)
w0=k_x*c/1
#print "w0",w0
a=wp_x_2/(1+p_y0_x)
#print('a',a)
'''
x_ref=ref_a(2800)
x_ref=np.array(x_ref)
print(x_ref)
######
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
plt.switch_backend('agg')
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
font = {'family' : 'Times New Roman',
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 15,  
        }  
color='rainbow'
index = 6

fig = plt.figure(1) 
#ax_cof = fig.add_subplot(121)

#fig = plt.figure(figsize=[9,6]) #定义figure，（1）中的1是什么
ax_cof = HostAxes(fig, [0, 0,0.9, 0.9])  #用[left, bottom, weight, height]的方式定义axes，0 <= l,b,w,h <= 1

#parasite addtional axes, share x
ax_temp = ParasiteAxes(ax_cof, sharex=ax_cof)
ax_load = ParasiteAxes(ax_cof, sharex=ax_cof)
ax_cp = ParasiteAxes(ax_cof, sharex=ax_cof)
ax_wear = ParasiteAxes(ax_cof, sharex=ax_cof)

#append axes
ax_cof.parasites.append(ax_temp)
ax_cof.parasites.append(ax_load)
ax_cof.parasites.append(ax_cp)
ax_cof.parasites.append(ax_wear)



#invisible right axis of ax_cof
ax_cof.axis['right'].set_visible(False)
ax_cof.axis['top'].set_visible(False)
ax_temp.axis['right'].set_visible(True)
ax_temp.axis['right'].major_ticklabels.set_visible(True)
ax_temp.axis['right'].label.set_visible(True)

#set label for axis
ax_cof.set_ylabel('k/${k_0}$')
ax_cof.set_xlabel('${/um}$')
ax_temp.set_ylabel('Ey')
ax_load.set_ylabel('refrective index')
#ax_cp.set_ylabel('k/${k_0}$')
#ax_wear.set_ylabel('Wear')

load_axisline = ax_load.get_grid_helper().new_fixed_axis
cp_axisline = ax_cp.get_grid_helper().new_fixed_axis
wear_axisline = ax_wear.get_grid_helper().new_fixed_axis

ax_load.axis['right2'] = load_axisline(loc='right', axes=ax_load, offset=(40,0))

#ax_cp.axis['right3'] = cp_axisline(loc='right', axes=ax_cp, offset=(80,0))
#ax_wear.axis['right4'] = wear_axisline(loc='right', axes=ax_wear, offset=(120,0))

fig.add_axes(ax_cof)

k0=2*3.14/const.lamada
k1=2*3.14/1e-6
fs=2*3.14/const.delta_x/k0
f,t,zxx=signal.stft(bz_y0,fs=2*3.14/const.delta_x,nperseg=200,noverlap=199)

plt.set_cmap(reg_cmap_transparent(color,create_alpha(lambda x:(np.exp(x/(256/5))-1.1))))
cmap = plt.get_cmap()

name='${e_k}$  '

curve_cof= ax_cof.pcolormesh(t*2*3.14/1e-6,f/k0,np.abs(zxx),cmap=cmap,shading='gouraud')


curve_temp= ax_temp.plot(np.linspace(0,60*10.6,1000),hy, label="Ey", color='red')
curve_load= ax_load.plot(np.linspace(0,60*10.6,1000),x_ref, label="ref", color='green')

ax_cp.set_ylim((0,2))



#ax_cof.set_ylim((-ne_y0.max(),ne_y0.max()))
ax_temp.set_ylim((-hy.max(),hy.max()))
ax_load.set_ylim((0,1))
ax_cof.set_ylim((0,2))


ax_cof.legend()

#轴名称，刻度值的颜色
#ax_cof.axis['left'].label.set_color(ax_cof.get_color())
ax_temp.axis['right'].label.set_color('red')
ax_load.axis['right2'].label.set_color('green')
#ax_cp.axis['right3'].label.set_color('pink')
#ax_wear.axis['right4'].label.set_color('blue')

ax_temp.axis['right'].major_ticks.set_color('red')
ax_load.axis['right2'].major_ticks.set_color('green')
#ax_cp.axis['right3'].major_ticks.set_color('pink')
#ax_wear.axis['right4'].major_ticks.set_color('blue')

ax_temp.axis['right'].major_ticklabels.set_color('red')
ax_load.axis['right2'].major_ticklabels.set_color('green')
#ax_cp.axis['right3'].major_ticklabels.set_color('pink')
#ax_wear.axis['right4'].major_ticklabels.set_color('blue')

ax_temp.axis['right'].line.set_color('red')
ax_load.axis['right2'].line.set_color('green')
#ax_cp.axis['right3'].line.set_color('pink')
#ax_wear.axis['right4'].line.set_color('blue')

position=fig.add_axes([0.05, 0.12, 0.28, 0.02])#位置[左,下,右,上]
cb=plt.colorbar(curve_cof,cax=position,orientation='horizontal')#方向
cb.ax.set_title('Wigner spectrogram', loc = 'left', fontdict=font)
cb.ax.tick_params(labelsize=15)
#plt.tight_layout()
plt.show()
fig.savefig(savedir,dpi=400,bbox_inches = 'tight')
