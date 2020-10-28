import numpy as np
import sdf
import constant as const
n=3000
sdfdir = const.sdfdir+str(n).zfill(4)+'.sdf'
savedir = const.txtdir + str(n)
data=sdf.read(sdfdir,dict=True)
bz=data['Electric Field/Ey'].data[:,:,int(188/2)]
ne=data['Derived/Number_Density/electron1'].data[:,:,int(188/2)]

np.savetxt(savedir+'_bz.txt',bz)
np.savetxt(savedir+'_ne.txt',ne)

