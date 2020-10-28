import numpy as np
import constant as const
a=np.zeros([4000,2])
for i in range(1,4000):
	a[i]=np.loadtxt("baktxt/"+const.data_name+str(i)+"E_2.txt")

b=a[:,1]/a[101,0]
print(b)
print(max(b))
