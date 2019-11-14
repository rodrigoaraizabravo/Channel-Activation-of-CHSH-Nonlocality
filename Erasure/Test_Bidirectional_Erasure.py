# -*- coding: utf-8 -*-
"""
@author: Rodrigo Araiza
"""
import Bidirectional_Erasure_Utils as BEU
import numpy as np
#from scipy.linalg import eigh

'''This scrip tests the bidirectional erasure functions'''

'''First a simple Identity channel test. We should get 2sqrt(2)'''
#K1=BEU.Erasure_Kraus(lam=0,where=1)
#K2=BEU.Erasure_Kraus(lam=0,where=2)
#K=[np.kron(k1,k2) for k2 in K2 for k1 in K1]

#x,states,meas=BEU.SeeSaw(100,20,K,stateFixed=False,track=True)
'''Test passed'''

'''Now, one channel is kept at p=0 while the other one is tempered with, 
one should still see max violation'''
#for p in np.linspace(0,1,11):
#    K1,K2=BEU.Erasure_Kraus(lam=0,where=1),BEU.Erasure_Kraus(lam=p,where=2)
#    K=[np.kron(k1,k2) for k2 in K2 for k1 in K1]
    
#    x,states,meas=BEU.SeeSaw(20,10,K,stateFixed=False,track=True)
#    print('\n')
'''Test passed'''
 
'''Monoticity test: one channel is set to 1/2 while the other one is free to vary'''
#CHSHVals=[]
#for p in np.linspace(0,1,11):
#    K1,K2=BEU.Erasure_Kraus(lam=1/2,where=1),BEU.Erasure_Kraus(lam=p,where=2)
#    K=[np.kron(k1,k2) for k2 in K2 for k1 in K1]
    
#    x,states,meas=BEU.SeeSaw(20,10,K,stateFixed=False,track=True)
#    CHSHVals.append(x)
#    print('\n') 
'''Test passed'''

'''Two erasure channels'''
#CHSHVals=np.zeros((21,21))
#ps=np.linspace(1/2,3/4,21)

#for i in range(len(ps)):
#    for j in range(len(ps)):
#        K1,K2=BEU.Erasure_Kraus(lam=ps[i],where=1),BEU.Erasure_Kraus(lam=ps[j],where=2)
#        K=[np.kron(k1,k2) for k2 in K2 for k1 in K1]
#        
#        x,states,meas=BEU.SeeSaw(30,50,K,stateFixed=False,track=True)
#        CHSHVals[i,j]=x
#        print('\n') 

'''Superactivation obeserved at p1=p2=1/2 2.0016'''


'''An erasure and a depolarization'''
#chshvals=[]
#for p in np.linspace(1-1/np.sqrt(2),1-1/np.sqrt(2),11):
#    K1,K2=BEU.Erasure_Kraus(lam=1/2,where=1),BEU.Depol_Kraus(lam=p,where=2)
#    K=[np.kron(k1,k2) for k2 in K2 for k1 in K1]
    
#    x,states,meas=BEU.SeeSaw(40,20,K,stateFixed=False,track=True)
#    chshvals.append(x)
#    print('\n') 


'''An erasure and an amplitude damping'''
#chshvals=[]
#for p in np.linspace(1/2,1/2,11):
#    K1,K2=BEU.Erasure_Kraus(lam=1/2,where=1),BEU.AmplitudeDamping_Kraus(lam=p,where=2)
#    K=[np.kron(k1,k2) for k2 in K2 for k1 in K1]
#    
#    x,states,meas=BEU.SeeSaw(40,20,K,stateFixed=False,track=True)
#    chshvals.append(x)
#    print('\n') 

pc=(3-np.sqrt(5))/2

'''An erasure and a lost channel'''
chshvals=[]
for p in np.linspace(0,pc,21):
    K1,K2=BEU.Erasure_Kraus(lam=1/2,where=1),BEU.Lost_Kraus(lam=p,where=2)
    K=[np.kron(k1,k2) for k2 in K2 for k1 in K1]
    x,states,meas=BEU.SeeSaw(40,40,K,stateFixed=False,track=True)
    chshvals.append(x)
    print('\n') 

