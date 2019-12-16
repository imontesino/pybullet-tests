import build.screwIK as IK
import numpy as np

x_com = 0
y_com = 0
z_com = 0.64

l0 = 0.1932
l1 = 0.305
l2 = 0.1625
l3 = 0.059742
l4 = 0.037508
l5 = 0.34692
l6 = 0.32992
l7 = 0.215
l8 = 0.090
l9 = 0.092
l10 = 0.330
l11 = 0.300
l12 = 0.123005
l13 = 0.146
l14 = 0.018
l15 = 0.026
l16 = 0.0175


# Right leg with respect to ground
xi1  = -np.array([1., 0., 0., 0., l12, 0.])
xi2  =  np.array([0., 1., 0., -l12, 0., 0.])
xi3  =  np.array([0., 1., 0., -(l11 + l12), 0., 0.])
xi4  =  np.array([0., 1., 0., -(l10 + l11 + l12), 0., 0.])
xi5  = -np.array([1., 0., 0., 0., (l10 + l11 + l12), (l16) ])
xi6  = -np.array([0., 0., 1., -l16, 0., 0.])


#Left leg with respect to to ground
xi7  =  np.array([0., 0., 1., l16, 0., 0.])
xi8  = -np.array([1., 0., 0., 0., (l10 + l11 + l12), -(l16) ])
xi9  =  np.array([0., 1., 0., -(l10 + l11 + l12), 0., 0.])
xi10 =  np.array([0., 1., 0., -(l11 + l12), 0., 0.])
xi11 =  np.array([0., 1., 0., -l12, 0., 0.])
xi12 = -np.array([1., 0., 0., 0., l12, 0.])

n=np.array([1, 0, 0])
o=np.array([0, 1, 0])
a=np.array([0, 0, 1])
point_p=np.array([x_com, y_com, z_com])-np.array([8.67/1000, -0.03/1000, 14.06/1000])+np.array([0, -l13, 0])-np.array(pos_r)

point_r=point_p

noap=Id(4)
noap[0:3,3]=point_p

Hst0=Id(4)
Hst0Inv=Id(4)

Hst0[0:3,3]=np.array([0, -l16, l9+l10+l11+l12])
Hst0Inv[0:3,3]=-np.array([0, -l16, l9+l10+l11+l12])

delta=((noap*Hst0Inv*cart2hom(f))[0:3,:]-k).norm()
theta3=symPK3(xi3,f,k,r,delta)[1]

print(IK.PK3(xi1,p,q,r, delta))
