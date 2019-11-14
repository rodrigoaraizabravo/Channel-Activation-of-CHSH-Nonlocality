# -*- coding: utf-8 -*-
"""
@author: Rodrigo Araiza
"""

'''This is a test on bidirectional protocol'''

import Superactivation_Utils as SU
import numpy as np
import matplotlib.pyplot as plt

def PTs(rho):
    rhoT=np.zeros(np.shape(rho),dtype=complex)
    rhoT[0,0],rhoT[0,1],rhoT[0,2],rhoT[0,3]=rho[0,0],rho[1,0],rho[0,2],rho[1,2]
    rhoT[1,0],rhoT[1,1],rhoT[1,2],rhoT[1,3]=rho[0,1],rho[2,2],rho[0,3],rho[1,3]
    rhoT[2,0],rhoT[2,1],rhoT[2,2],rhoT[2,3]=rho[2,0],rho[3,0],rho[2,2],rho[3,2]
    rhoT[3,0],rhoT[3,1],rhoT[3,2],rhoT[3,3]=rho[2,1],rho[3,1],rho[2,3],rho[3,3]
    return rhoT

'''Two bit flip channels'''
'''
psamples=np.linspace(0,1,11)
for p in psamples:
    K1,K2=SU.BitFlip_Kraus(lam=p,where=1),SU.BitFlip_Kraus(lam=p,where=2)
    print(p)
    xtop=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)
print(xtop)'''

'''Two phase flip channels'''
'''
psamples=np.linspace(0,1,21)
for p in psamples:
    K1,K2=SU.PhaseFlip_Kraus(lam=p,where=1),SU.PhaseFlip_Kraus(lam=p,where=2)
    xtop=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)'''

'''A bit and a phase flip channels'''
'''
psamples=np.linspace(1/2,1/2,1)
for p in psamples:
    K1,K2=SU.BitFlip_Kraus(lam=p,where=1),SU.PhaseFlip_Kraus(lam=p,where=2)
    xtop=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)'''

'''Two depolarization channels'''
'''
p1samples=np.linspace(1-1/np.sqrt(2),1,11)
p2samples=np.linspace(0,1,11)
for p1 in p1samples:
    for p2 in p2samples:
        K1,K2=SU.Depol_Kraus(lam=p1,where=1),SU.Depol_Kraus(lam=p2,where=2)
        xtop=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)'''

'''One depolarization and one bit flip channels'''
'''
p1samples=np.linspace(1-1/np.sqrt(2),1,11)
p2samples=np.linspace(0,1/2,11)
for p1 in p1samples:
    for p2 in p2samples:
        K1,K2=SU.Depol_Kraus(lam=p1,where=1),SU.BitFlip_Kraus(lam=p2,where=2)
        xtop=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)'''

'''One depolarization and one phase flip channels'''
'''
p1samples=np.linspace(1-1/np.sqrt(2),1,11)
p2samples=np.linspace(1/2,1/2,1)
for p1 in p1samples:
    for p2 in p2samples:
        K1,K2=SU.Depol_Kraus(lam=p1,where=1),SU.PhaseFlip_Kraus(lam=p2,where=2)
        xtop=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)'''

'''Two amplitude damping channels'''
'''
p1samples=np.linspace(1/2,1/2,1)
p2samples=np.linspace(1/2,1/2,1)
CHSHVals=np.zeros((len(p1samples),len(p2samples)))
for i in range(len(p1samples)):
    for j in range(len(p2samples)):
        p1,p2=p1samples[i],p2samples[j]
        K1,K2=SU.AmplitudeDamping_Kraus(lam=p1,where=1),SU.AmplitudeDamping_Kraus(lam=p2,where=2)
        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)
        CHSHVals[i,j]=xtop
'''
'''
plt.figure()
color_map = plt.imshow(CHSHVals[10:,10:], extent=[0.5,0.7,0.7,0.5])
color_map.set_cmap("cool")
plt.colorbar()
plt.xlabel(r'$p_1$',fontsize=16)
plt.ylabel(r'$p_2$',fontsize=16)
plt.title('CHSH superactivation of amp. damp. channels. \n Max=2.0119',fontsize=18)
plt.show()
'''
'''Superactivation observed at p1=[0.5,0.55] and p2=[1/2,0.65]'''

'''One amplitude damping channel and one bit flip'''
'''
p1samples=np.linspace(1/2,1,11)
p2samples=np.linspace(0,1/2,11)
for p1 in p1samples:
    for p2 in p2samples:
        K1,K2=SU.AmplitudeDamping_Kraus(lam=p1,where=1),SU.BitFlip_Kraus(lam=p2,where=2)
        xtop=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)'''

'''One amplitude damping channel and one phase flip'''
'''
p1samples=np.linspace(1/2,1,11)
p2samples=np.linspace(0,1/2,11)
for p1 in p1samples:
    for p2 in p2samples:
        K1,K2=SU.AmplitudeDamping_Kraus(lam=p1,where=1),SU.PhaseFlip_Kraus(lam=p2,where=2)
        xtop=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)'''
'''No superactivation'''

'''One amplitude damping channel and one depolarization channel'''
'''
p1samples=np.linspace(1/2,1/2,1)
p2samples=np.linspace(1-1/np.sqrt(2),1-1/np.sqrt(2),1)
for p1 in p1samples:
    for p2 in p2samples:
        K1,K2=SU.AmplitudeDamping_Kraus(lam=p1,where=1),SU.Depol_Kraus(lam=p2,where=2)
        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,iters=40,track=True)
'''

'''Two loss channels'''
pc=(3-np.sqrt(5))/2
'''
p1samples=np.linspace(pc,pc,6)
p2samples=np.linspace(pc,pc,1)
for p1 in p1samples:
    for p2 in p2samples:
        K1,K2=SU.Loss_Kraus(lam=p1,where=1),SU.Loss_Kraus(lam=p2,where=2)
        xtop=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=True)
'''

'''Loss channel and bit flip'''
'''
p1samples=np.linspace(pc,pc,11)
p2samples=np.linspace(1/2,1/2,1)
xtop=0
for p1 in p1samples:
    for p2 in p2samples:
        K1,K2=SU.Loss_Kraus(lam=p1,where=1),SU.BitFlip_Kraus(lam=p2,where=2)
        x,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=True)
        if x>xtop:xtop=x
'''

'''Loss channel and depolarization channel'''
'''
p1samples=np.linspace(1/2,1/2,50)
p2samples=np.linspace(1-1/np.sqrt(2),1-1/np.sqrt(2),1)
xtop=0
for p1 in p1samples:
    for p2 in p2samples:
        K1,K2=SU.Loss_Kraus(lam=p1,where=1),SU.Depol_Kraus(lam=p2,where=2)
        x,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1,pJump=1,iters=40,track=True)
        if x>xtop:xtop=x
'''

'''Loss channel and amplitude damping channel'''

p1samples=np.linspace(1/2,1/2,50)
p2samples=np.linspace(1/2,1/2,1)
xtop=0
for p1 in p1samples:
    for p2 in p2samples:
        K1,K2=SU.Lost_Kraus(lam=p1,where=1),SU.AmplitudeDamping_Kraus(lam=p2,where=2)
        x,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1,pJump=1,iters=40,track=True)
        if x>xtop:xtop=x
  
'''Activation of 2.00212 at pc'''
'''Activation of 2.00023'''

