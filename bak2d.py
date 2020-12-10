import numpy as np
import sdf
import constant as const
n=1600
sdfdir = const.sdfdir+str(n).zfill(4)+'.sdf'
savedir = const.txtdir + str(n)
data=sdf.read(sdfdir,dict=True)
Ex=data["Electric Field/Ex"].data[:,:,int(188/2)]
Ey=data['Electric Field/Ey'].data[:,:,int(188/2)]
Bz=data['Magnetic Field/Bz'].data[:,:,int(188/2)]
ne=data['Derived/Number_Density/electron1'].data[:,:,int(188/2)]
np.savetxt(savedir+'_ex.txt',Ex)
np.savetxt(savedir+'_ey.txt',Ey)
np.savetxt(savedir+'_bz.txt',Bz)
np.savetxt(savedir+'_ne.txt',ne)

