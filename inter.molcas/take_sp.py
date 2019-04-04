import h5py,sys
import numpy as np

file = sys.argv[1]
with h5py.File( file ) as f:
# for key in f:
#  print(key)
# print(f['CENTER_CHARGES'][:])
# print(f['SOS_ENERGIES'])
 E=f['SFS_ENERGIES'][:]
 dip=f['SFS_EDIPMOM'][:]

npot=len(E)

Hr=np.zeros( (npot,npot) )
for i in range(npot):
 Hr[i,i]=E[i]

np.savetxt("H0.dat",Hr)

#### CREATING DIPOLE MATRIX
dip1=dip[0,:,:]
dip2=dip[1,:,:]
dip3=dip[2,:,:]

np.savetxt("dipole1.dat",dip1)
np.savetxt("dipole2.dat",dip2)
np.savetxt("dipole3.dat",dip3)
