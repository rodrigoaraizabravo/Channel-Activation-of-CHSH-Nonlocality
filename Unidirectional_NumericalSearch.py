# -*- coding: utf-8 -*-
"""
@author: Rodrigo Araiza
"""
'''This is a test on uni-directional protocol'''

import Superactivation_Utils as SU
import numpy as np

'''Two bit flip channels'''
psamples=np.linspace(0,0.2,2)
for p in psamples:
    K1,K2=SU.BitFlip_Kraus(lam=p,where=1),SU.BitFlip_Kraus(lam=p,where=1)
    print(p)
    xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)
print(xtop)
'''No superactivation'''

'''Two phase flip channels'''
psamples=np.linspace(0,1,21)[:-1]
for p in psamples:
    K1,K2=SU.PhaseFlip_Kraus(lam=p,where=1),SU.PhaseFlip_Kraus(lam=p,where=1)
    xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)
'''No superactivation'''

'''A bit and a phase flip channels'''
psamples=np.linspace(0,0.50,11)
for p in psamples:
    K1,K2=SU.BitFlip_Kraus(lam=p,where=1),SU.PhaseFlip_Kraus(lam=p,where=1)
    xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)
'''No superactivation'''

'''Two depolarization channels'''
#p1samples=np.linspace(1-1/np.sqrt(2),1,11)
#p2samples=np.linspace(0,1-1/np.sqrt(2),11)
#for p1 in p1samples:
#    for p2 in p2samples:
#        K1,K2=SU.Depol_Kraus(lam=p1,where=1),SU.Depol_Kraus(lam=p2,where=1)
#        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1/2,pJump=1/2,track=False)
'''No superactivation'''

'''One depolarization and one bit flip channels'''
#p1samples=np.linspace(1-1/np.sqrt(2),1,11)
#p2samples=np.linspace(0,1/2,11)
#for p1 in p1samples:
#    for p2 in p2samples:
#        K1,K2=SU.Depol_Kraus(lam=p1,where=1),SU.BitFlip_Kraus(lam=p2,where=1)
#        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)
'''No superactivation'''

'''One depolarization and one phase flip channels'''
#p1samples=np.linspace(1-1/np.sqrt(2),1,11)
#p2samples=np.linspace(1/2,1/2,1)
#for p1 in p1samples:
#    for p2 in p2samples:
#        K1,K2=SU.Depol_Kraus(lam=p1,where=1),SU.PhaseFlip_Kraus(lam=p2,where=1)
#        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1e-2,pJump=1/2,track=False)
'''No superactivation'''

'''Two amplitude damping channels'''
#p1samples=np.linspace(1/2,1,11)
#p2samples=np.linspace(1/2,1,11)
#for p1 in p1samples:
#    for p2 in p2samples:
#        K1,K2=SU.AmplitudeDamping_Kraus(lam=p1,where=1),SU.AmplitudeDamping_Kraus(lam=p2,where=1)
#        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1/2,maxtimes=30,pJump=1/2,track=False)
'''No superactivation observed'''

'''One amplitude damping channel and one bit flip'''
#p1samples=np.linspace(1/2,1,11)
#p2samples=np.linspace(1/2,1/2,1)
#for p1 in p1samples:
#    for p2 in p2samples:
#        K1,K2=SU.AmplitudeDamping_Kraus(lam=p1,where=1),SU.BitFlip_Kraus(lam=p2,where=1)
#        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1/2,maxtimes=30,pJump=1/2,track=False)
'''No superactivation'''

'''One amplitude damping channel and one phase flip'''
#p1samples=np.linspace(1/2,1,11)
#p2samples=np.linspace(1/2,1/2,1)
#for p1 in p1samples:
#    for p2 in p2samples:
#        K1,K2=SU.AmplitudeDamping_Kraus(lam=p1,where=1),SU.PhaseFlip_Kraus(lam=p2,where=1)
#        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1/2,maxtimes=30,pJump=1/2,track=False)
#'''No superactivation'''

'''One amplitude damping channel and one depolarization flip'''
#p1samples=np.linspace(1/2,1,11)
#p2samples=np.linspace(1-1/np.sqrt(2),1,11)
#for p1 in p1samples:
#    for p2 in p2samples:
#        K1,K2=SU.AmplitudeDamping_Kraus(lam=p1,where=1),SU.Depol_Kraus(lam=p2,where=1)
#        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1/2,maxtimes=30,pJump=1/2,track=False)
'''No superactivation'''
pc=0.34130739005137906
'''Two loss channels'''
#p1samples=np.linspace(pc,pc,11)
#p2samples=np.linspace(pc,pc,1)
#for p1 in p1samples:
#    for p2 in p2samples:
#        K1,K2=SU.Loss_Kraus(lam=p1,where=1),SU.Lost_Kraus(lam=p2,where=1)
#        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1/2,maxtimes=30,pJump=1/2,track=True)
'''Activation of 2.02618 '''

'''A loss channel and a bit flip channel'''
#p1samples=np.linspace(pc,pc,11)
#p2samples=np.linspace(1/2,1/2,1)
#for p1 in p1samples:
#    for p2 in p2samples:
#        K1,K2=SU.Loss_Kraus(lam=p1,where=1),SU.BitFlip_Kraus(lam=p2,where=1)
#        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1/2,maxtimes=30,pJump=1/2,track=True)
'''Also activation of 2.026818'''

'''A loss channel and a depolarization channel'''
#p1samples=np.linspace(pc,pc,11)
#p2samples=np.linspace(1-1/np.sqrt(2),1-1/np.sqrt(2),1)
#for p1 in p1samples:
#    for p2 in p2samples:
#        K1,K2=SU.Loss_Kraus(lam=p1,where=1),SU.Depol_Kraus(lam=p2,where=1)
#        xtop,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,eps=1/2,maxtimes=30,pJump=1/2,track=True)
'''No superactivation'''
'''Also activation of 2.02681'''

'''A loss channel and an amplitude damping'''
#xtop=0
#p1samples=np.linspace(pc,pc,1)
#p2samples=np.linspace(1/2,1,11)
#for p1 in p1samples:
#    for p2 in p2samples:
#       K1,K2=SU.Loss_Kraus(lam=p1,where=1),SU.AmplitudeDamping_Kraus(lam=p2,where=1)
#       x,state=SU.ChannelsSeeSaw(1e-5,20,K1,K2,iters=40,eps=1,maxtimes=2,pJump=0,track=True)
#       if x>xtop: xtop=x
'''Also activation of 2.02681'''