# -*- coding: utf-8 -*-
"""
@author: Rodrigo Araiza
"""

import numpy as np
import Erasure_Utils as EU
from scipy.linalg import eigh


'''A amplitude damping channel and a depolarization channel full state input:
    WE SHOULD see superactivation'''
'''
for p in np.linspace(1/2,1/2,1):
    K1,K2=EU.AmplitudeDamping_Kraus(lam=p),EU.Depol_Kraus(lam=1-1/np.sqrt(2))
    K=[np.kron(K1[i],K2[j]) for i in range(len(K1)) for j in range(len(K2))]
    B=np.kron(np.kron(EU.QutritMeasRand(),EU.QubitMeasRand()),np.kron(EU.QutritMeasRand(),EU.QubitMeasRand()))
    ev,evecs=eigh(B)
    rho=1/2*(np.outer(evecs[:,-1],np.conj(evecs[:,-1]))+np.outer(evecs[:,-2],np.conj(evecs[:,-2])))
    x,state=EU.UntamperedSeeSawFullState(40,20,rho,K,stateFixed=False,track=True)
    print('\n')
'''

'''Two erasure channels''' 
'''
for p in np.linspace(0.5,0.5,11):
    K1,K2=EU.Erasure_Kraus(lam=p),EU.Erasure_Kraus(lam=1/2)
    K=[np.kron(K1[i],K2[j]) for j in range(len(K2)) for i in range(len(K1))]
    B=np.kron(np.kron(EU.QutritMeasRand(),EU.QubitMeasRand()),np.kron(EU.QutritMeasRand(),EU.QubitMeasRand()))
    ev,evecs=eigh(B)
    rho=1/2*(np.outer(evecs[:,-1],np.conj(evecs[:,-1]))+np.outer(evecs[:,-2],np.conj(evecs[:,-2])))
    x,state=EU.UntamperedSeeSawFullState(100,10,rho,K,stateFixed=False,track=True)
    if x<2.0: break
    #x,state=EU.SeeSawFullState(1e-5,100,20,rho,K,eps=1/2,pJump=1/2,stateFixed=False,track=True)
    print('\n')
'''
'''No superactivation observed'''

'''One erasure channel and one amplitude damping channel''' 
'''
for p in np.linspace(0,0.5,11):
    K1,K2=EU.Erasure_Kraus(lam=p),EU.AmplitudeDamping_Kraus(lam=1/2)
    K=[np.kron(K1[i],K2[j]) for j in range(len(K2)) for i in range(len(K1))]
    B=np.kron(np.kron(EU.QutritMeasRand(),EU.QubitMeasRand()),np.kron(EU.QutritMeasRand(),EU.QubitMeasRand()))
    ev,evecs=eigh(B)
    rho=1/2*(np.outer(evecs[:,-1],np.conj(evecs[:,-1]))+np.outer(evecs[:,-2],np.conj(evecs[:,-2])))
    x,state=EU.UntamperedSeeSawFullState(100,10,rho,K,stateFixed=False,track=True)
    if x<2.0: break
    print('\n')
    '''
    
'''One erasure channel and one depolarizing channel''' 
'''
for p in np.linspace(0.5,0.5,11):
    xtop=2.0
    K1,K2=EU.Erasure_Kraus(lam=p),EU.Depol_Kraus(lam=1-1/np.sqrt(2))
    K=[np.kron(K1[i],K2[j]) for j in range(len(K2)) for i in range(len(K1))]
    B=np.kron(np.kron(EU.QutritMeasRand(),EU.QubitMeasRand()),np.kron(EU.QutritMeasRand(),EU.QubitMeasRand()))
    ev,evecs=eigh(B)
    rho=1/2*(np.outer(evecs[:,-1],np.conj(evecs[:,-1]))+np.outer(evecs[:,-2],np.conj(evecs[:,-2])))
    x,state=EU.UntamperedSeeSawFullState(100,20,rho,K,stateFixed=False,track=True)
    if x>xtop: xtop=x
    print('\n')
'''

'''Loss channel and erasure'''
pc=(3-np.sqrt(5))/2

for p in np.linspace(pc,pc,11):
    xtop=2.0
    K1,K2=EU.Loss_Kraus(lam=p),EU.Erasure_Kraus(lam=1/2)
    K=[np.kron(K1[i],K2[j]) for j in range(len(K2)) for i in range(len(K1))]
    B=np.kron(np.kron(EU.QutritMeasRand(),EU.QubitMeasRand()),np.kron(EU.QutritMeasRand(),EU.QubitMeasRand()))
    ev,evecs=eigh(B)
    rho=1/2*(np.outer(evecs[:,-1],np.conj(evecs[:,-1]))+np.outer(evecs[:,-2],np.conj(evecs[:,-2])))
    x,state=EU.UntamperedSeeSawFullState(100,20,rho,K,stateFixed=False,track=True)
    if x>xtop: xtop=x
    print('\n')