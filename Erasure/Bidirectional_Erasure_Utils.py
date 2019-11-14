# -*- coding: utf-8 -*-
"""
@author: Rodrigo Araiza
"""

'''This is a library of functions for the erasure channel in the bi-directional
protocol'''

import numpy as np
from scipy.linalg import eigh
from Erasure_Utils import ChannelState, ConstBp, QubitMeasRand, QutritMeasRand, Block

'''We imagine the erasure channel acting on the first qubit thus we have a qutrit-qubit 
system'''
'''Hardcoding some important states and matrices'''
GHZ=np.array([0,1,-1,0,0,0])/np.sqrt(2)
rhoGHZ=np.outer(GHZ,GHZ)
sigma2z,sigma3z=np.array([[1,0],[0,-1]]),np.array([[1,0,0],[0,-1,0],[0,0,1]])
sigma2x,sigma3x=np.array([[0,1],[1,0]]),np.array([[0,1,0],[1,0,0],[0,0,1]])
sigma2y,sigma3y=np.array([[0,-1j],[1j,0]]),np.array([[0,-1j,0],[1j,0,0],[0,0,1]])
Id9,Id6,Id4,Id3,Id2=np.identity(9),np.identity(6),np.identity(4),np.identity(3),np.identity(2)

def Erasure_Kraus(lam=1/2,where=1):
    '''Kraus ops for erasure channel'''
    if where==1:
        K =[np.sqrt(1-lam)*np.identity(6)]
        K+=[np.sqrt(lam)*np.kron(np.outer(np.array([0,0,1]),np.array([1,0,0])),Id2)]
        K+=[np.sqrt(lam)*np.kron(np.outer(np.array([0,0,1]),np.array([0,1,0])),Id2)]
        K+=[np.sqrt(lam)*np.kron(np.outer(np.array([0,0,1]),np.array([0,0,1])),Id2)]
    
    elif where==2:
        K =[np.sqrt(1-lam)*np.identity(6)]
        K+=[np.sqrt(lam)*np.kron(Id2,np.outer(np.array([0,0,1]),np.array([1,0,0])))]
        K+=[np.sqrt(lam)*np.kron(Id2,np.outer(np.array([0,0,1]),np.array([0,1,0])))]
        K+=[np.sqrt(lam)*np.kron(Id2,np.outer(np.array([0,0,1]),np.array([0,0,1])))]
    return K

def Loss_Kraus(lam=1/2,where=1):
    '''Kraus ops for loss channel'''
    if where==1:
        K =[np.sqrt(1-lam)*np.identity(6)]
        K+=[np.sqrt(lam)*np.kron(np.outer(np.array([1,0,0]),np.array([1,0,0])),Id2)]
        K+=[np.sqrt(lam)*np.kron(np.outer(np.array([1,0,0]),np.array([0,1,0])),Id2)]
        K+=[np.sqrt(lam)*np.kron(np.outer(np.array([1,0,0]),np.array([0,0,1])),Id2)]
        
    elif where==2:
        K =[np.sqrt(1-lam)*np.identity(6)]
        K+=[np.sqrt(lam)*np.kron(Id2,np.outer(np.array([1,0,0]),np.array([1,0,0])))]
        K+=[np.sqrt(lam)*np.kron(Id2,np.outer(np.array([1,0,0]),np.array([0,1,0])))]
        K+=[np.sqrt(lam)*np.kron(Id2,np.outer(np.array([1,0,0]),np.array([0,0,1])))]
    return K  

def AmplitudeDamping_Kraus(lam=1/2,where=1):
    '''Kraus ops for amplitude damping channel'''
    if where==1:
        K =[np.kron(np.array([[1,0,0],[0,np.sqrt(1-lam),0],[0,0,1]]),Id2)]
        K+=[np.kron(np.array([[0,np.sqrt(lam),0],[0,0,0],[0,0,0]]),Id2)]
    elif where==2:
        K =[np.kron(Id2,np.array([[1,0,0],[0,np.sqrt(1-lam),0],[0,0,1]]))]
        K+=[np.kron(Id2,np.array([[0,np.sqrt(lam),0],[0,0,0],[0,0,0]]))]
    return K

def Depol_Kraus(lam=1-1/np.sqrt(2),where=1):
    '''Kraus operators for the depolarization channel'''
    if where==1:
        K =[np.sqrt(1-3/4*lam)*np.identity(6)]
        K+=[np.sqrt(lam/4)*np.kron(sigma3x,np.identity(2))]
        K+=[np.sqrt(lam/4)*np.kron(sigma3y,np.identity(2))]
        K+=[np.sqrt(lam/4)*np.kron(sigma3z,np.identity(2))]
    elif where==2:
        K =[np.sqrt(1-3/4*lam)*np.identity(6)]
        K+=[np.sqrt(lam/4)*np.kron(Id2,sigma3x)]
        K+=[np.sqrt(lam/4)*np.kron(Id2,sigma3y)]
        K+=[np.sqrt(lam/4)*np.kron(Id2,sigma3z)]
    return K

pF=np.argsort([0,1,4,2,3,5])
pB=np.argsort([0,1,3,4,2,5])

def SwitchNsBasis(N,mode='F'):
    '''from a matrix, it makes anotherone acting on a larger space'''
    if mode=='F': 
        s=N[:,pF]
        return s[pF,:]
    elif mode=='B': 
        s=N[pB,:]
        return s[:,pB]

'''Trace functions'''
def TrAAp(m):
    rho1=np.trace(m.reshape([3,2,2,3,3,2,2,3]),axis1=0,axis2=4)
    return np.trace(rho1,axis1=1,axis2=4).reshape((6,6))

def TrBBp(m):
    rho1=np.trace(m.reshape([3,2,2,3,3,2,2,3]),axis1=1,axis2=5)
    return np.trace(rho1,axis1=2,axis2=5).reshape((6,6))

def TrApBp(m):
    rho1=np.trace(m.reshape([3,2,2,3,3,2,2,3]),axis1=2,axis2=6)
    return np.trace(rho1,axis1=2,axis2=5).reshape((6,6))

def TrAB(m):
    rho1=np.trace(m.reshape([3,2,2,3,3,2,2,3]),axis1=0,axis2=4)
    return np.trace(rho1,axis1=0,axis2=3).reshape((6,6))

'''Functions to find optimal measurement settings'''
def FindMs(N1,N2,state):
    F1=TrBBp((N1+N2)@state) #6x6
    m1=np.zeros((6,6),dtype=complex)
    ev11,evecs11=eigh(F1[0:4,0:4])
    ev12,evecs12=eigh(F1[4:6,4:6])
    m1[0:4,0:4],m1[4:6,4:6]=Block(ev11,evecs11),Block(ev12,evecs12)
    
    F2=TrBBp((N1-N2)@state) #6x6
    ev21,evecs21=eigh(F2[0:4,0:4])
    ev22,evecs22=eigh(F2[4:6,4:6])
    m2=np.zeros((6,6),dtype=complex)
    m2[0:4,0:4],m2[4:6,4:6]=Block(ev21,evecs21),Block(ev22,evecs22)

    return m1,m2

def FindNs(M1,M2,state):
    F1=SwitchNsBasis(TrAAp((M1+M2)@state),mode='F')
    n1=np.zeros((6,6),dtype=complex)
    ev11,evecs11=eigh(F1[0:4,0:4])
    ev12,evecs12=eigh(F1[4:6,4:6])
    n1[0:4,0:4],n1[4:6,4:6]=Block(ev11,evecs11),Block(ev12,evecs12)
    
    F2=SwitchNsBasis(TrAAp((M1-M2)@state),mode='F')
    ev21,evecs21=eigh(F2[0:4,0:4])
    ev22,evecs22=eigh(F2[4:6,4:6])
    n2=np.zeros((6,6),dtype=complex)
    n2[0:4,0:4],n2[4:6,4:6]=Block(ev21,evecs21),Block(ev22,evecs22)
    
    return SwitchNsBasis(n1,mode='B'),SwitchNsBasis(n2,mode='B')
'''Functions to find optimal states from a bell operator'''
def Findrho1(Bp,rho2):
    ev,evec = eigh(TrApBp(np.kron(Id6,rho2)@Bp)[0:4,0:4])
    rho1=np.zeros((6,6),dtype=complex)
    rho1[0:4,0:4] =np.outer(evec[:,-1],np.conj(evec[:,-1]))
    return rho1

def Findrho2(Bp,rho1):
    ev,evec = eigh(SwitchNsBasis(TrAB(np.kron(rho1,Id6)@Bp),mode='F')[0:4,0:4])
    rho2=np.zeros((6,6),dtype=complex)
    rho2[0:4,0:4] =np.outer(evec[:,-1],np.conj(evec[:,-1]))
    return SwitchNsBasis(rho2,mode='B')

keys=[a+b+ap+bp for a in ['0','1','2'] for b in ['0','1'] for ap in ['0','1'] for bp in ['0','1','2']]
keysM=['00','01','10','11','20','21']
keysN=['00','01','02','10','11','12']

def Pad(Key,m):
    '''Given a matrix acting locally makes a global matrix of the form I\otimesM\otimes I
    '''
    if Key=='M':
        O=np.zeros((36,36),dtype=complex)
        for i in range(36):
            for j in range(36):
                srow,scol=keys[i],keys[j]
                if srow[1]==scol[1] and srow[3]==scol[3]:
                    O[i,j]=m[keysM.index(srow[0]+srow[2]),keysM.index(scol[0]+scol[2])]
    elif Key=='N':
        O=np.zeros((36,36),dtype=complex)
        for i in range(36):
            for j in range(36):
                srow,scol=keys[i],keys[j]
                if srow[0]==scol[0] and srow[2]==scol[2]:
                    O[i,j]=m[keysN.index(srow[1]+srow[3]),keysN.index(scol[1]+scol[3])]
    return O

def RandState():
    '''Produces a random state'''
    B=np.kron(QutritMeasRand(),QubitMeasRand())
    ev1,evecs1=eigh(B)
    return np.outer(evecs1[:,-1],np.conj(evecs1[:,-1]))

def SeeSaw(iters,epochs,K,stateFixed=False,track=True):
    '''aldorithm'''
    xtop=0
    for epoch in range(epochs):
        M1=Pad('M',np.kron(QutritMeasRand(),QubitMeasRand()))
        M2=Pad('M',np.kron(QutritMeasRand(),QubitMeasRand()))
        rho1,rho2=RandState(),RandState()
        state=ChannelState(np.kron(rho1,rho2),K)
                
        for i in range(iters):
            n1,n2=FindNs(M1,M2,state)
            N1,N2=Pad('N',n1),Pad('N',n2)
            if stateFixed==False:
                rho1=Findrho1(ConstBp(M1@(N1+N2)+M2@(N1-N2),K),rho2)
                state=ChannelState(np.kron(rho1,rho2),K)
            m1,m2=FindMs(N1,N2,state)
            M1,M2=Pad('M',m1),Pad('M',m2)
            if stateFixed==False:
                rho2=Findrho2(ConstBp(M1@(N1+N2)+M2@(N1-N2),K),rho1)
                state=ChannelState(np.kron(rho1,rho2),K)
                    
        x1=np.abs(np.trace((M1@(N1+N2)+M2@(N1-N2))@state))

        if x1>xtop:
            xtop=x1
            fstate=[rho1,rho2]
            measurements=[m1,m2,n1,n2]
        if track==True: print(xtop)
    print(xtop)
    return xtop,fstate,measurements

