# -*- coding: utf-8 -*-
"""
@author: Rodrigo Araiza
"""

'''This is a library of functions having to do with the erasure channel'''
import numpy as np
from scipy.linalg import eigh
from scipy.stats import unitary_group
import sys
sys.path.insert(0, '/Users/oscar/OneDrive/Desktop/Superactivation/Channels')
from Superactivation_Utils import sigmaz


'''We imagine the erasure channel acting on the first qubit thus we have a qutrit-qubit 
system'''
GHZ=np.array([0,1,-1,0,0,0])/np.sqrt(2)
rhoGHZ=np.outer(GHZ,GHZ)
sigma2z,sigma3z=np.array([[1,0],[0,-1]]),np.array([[1,0,0],[0,-1,0],[0,0,1]])
sigma2x,sigma3x=np.array([[0,1],[1,0]]),np.array([[0,1,0],[1,0,0],[0,0,1]])
sigma2y,sigma3y=np.array([[0,-1j],[1j,0]]),np.array([[0,-1j,0],[1j,0,0],[0,0,1]])

def Erasure_Kraus(lam=1/2,where=1):
    if where!=1:raise ValueError('Erasure channel can only act on the first qutrit')
    K =[np.sqrt(1-lam)*np.identity(6)]
    K+=[np.sqrt(lam)*np.kron(np.outer(np.array([0,0,1]),np.array([1,0,0])),np.identity(2))]
    K+=[np.sqrt(lam)*np.kron(np.outer(np.array([0,0,1]),np.array([0,1,0])),np.identity(2))]
    K+=[np.sqrt(lam)*np.kron(np.outer(np.array([0,0,1]),np.array([0,0,1])),np.identity(2))]
    return K

def Loss_Kraus(lam=0.5,where=1,d=2):
    K =[Id6*np.sqrt(1-lam)]
    if where==1:
        K+=[np.sqrt(lam)*np.kron(np.array([[1,0,0],[0,0,0],[0,0,0]]),Id2)]
        K+=[np.sqrt(lam)*np.kron(np.array([[0,1,0],[0,0,0],[0,0,0]]),Id2)]
    else: raise ValueError('Lost channel can only act on first qutrit')
    return K

def AmplitudeDamping_Kraus(lam=1/2,where=1):
    if where!=1:raise ValueError('AmpDamp channel can only act on the first qutrit')
    K =[np.kron(np.array([[1,0,0],[0,np.sqrt(1-lam),0],[0,0,1]]),Id2)]
    K+=[np.kron(np.array([[0,np.sqrt(lam),0],[0,0,0],[0,0,0]]),Id2)]
    return K

def Depol_Kraus(lam=1-1/np.sqrt(2),where=1):
    if where!=1:raise ValueError('Erasure channel can only act on the first qutrit')
    K =[np.sqrt(1-3/4*lam)*np.identity(6)]
    K+=[np.sqrt(lam/4)*np.kron(sigma3x,np.identity(2))]
    K+=[np.sqrt(lam/4)*np.kron(sigma3y,np.identity(2))]
    K+=[np.sqrt(lam/4)*np.kron(sigma3z,np.identity(2))]
    return K
    
def ChannelState(rho,K):
    O=np.zeros(np.shape(rho),dtype=complex)
    for k in K:
        O+=k@rho@np.conj(k.T)
    return O

def TrA(m):
    '''Partial trace over the qutrit system leaving it a qubit system'''
    return np.trace(m.reshape([3,2,3,2]),axis1=0,axis2=2).reshape((2,2))

def TrB(m):
    '''Partial trace over the qubit system leaving it a qutrit system'''
    return np.trace(m.reshape([3,2,3,2]),axis1=1,axis2=3).reshape((3,3))

def QubitMeasRand():
    U=unitary_group.rvs(2)
    return U@sigmaz@np.conj(U.T)

def QutritMeasRand():
    M=np.zeros((3,3),dtype=complex)
    M[0:2,0:2]=QubitMeasRand()
    M[2,2]=2*np.random.binomial(1,0.5)-1
    return M

def FindMsErasure(N1,N2,rho):
    F1,F2=TrB(np.kron(np.identity(3),N1+N2)@rho),TrB(np.kron(np.identity(3),N1-N2)@rho)
    ev1,evec1=eigh(F1)
    ev2,evec2=eigh(F2)
    return evec1@np.diag(np.sign(ev1))@np.conj(evec1.T), evec2@np.diag(np.sign(ev2))@np.conj(evec2.T)

def FindNsErasure(M1,M2,rho):
    F1,F2=TrA(np.kron(M1+M2,np.identity(2))@rho),TrA(np.kron(M1-M2,np.identity(2))@rho)
    ev1,evec1=eigh(F1)
    ev2,evec2=eigh(F2)
    return evec1@np.diag(np.sign(ev1))@np.conj(evec1.T), evec2@np.diag(np.sign(ev2))@np.conj(evec2.T)

def MaxCHSHVal(K,iters=20,rounds=20,track=True,stateFixed=True):
    xtop=0.0
    for epoch in range(rounds):
        M1,M2,N1,N2=QutritMeasRand(),QutritMeasRand(),QubitMeasRand(),QubitMeasRand()
        rho=np.zeros((6,6),dtype=complex)
        state=ChannelState(rho,K)
        ev,evecs=eigh((np.kron(M1,N1+N2)+np.kron(M2,N1-N2))[0:4,0:4])
        rho[0:4,0:4]=np.outer(evecs[:,-1],np.conj(evecs[:,-1]))
        state=ChannelState(rho,K)
        
        for iter in range(iters):
            M1,M2=FindMsErasure(N1,N2,state)
            if stateFixed==False:
                ev,evecs=eigh((np.kron(M1,N1+N2)+np.kron(M2,N1-N2))[0:4,0:4])
                rho[0:4,0:4]=np.outer(evecs[:,-1],np.conj(evecs[:,-1]))
                state=ChannelState(rho,K)
            N1,N2=FindNsErasure(M1,M2,state)
            if stateFixed==False:
                ev,evecs=eigh((np.kron(M1,N1+N2)+np.kron(M2,N1-N2))[0:4,0:4])
                rho[0:4,0:4]=np.outer(evecs[:,-1],np.conj(evecs[:,-1]))
                state=ChannelState(rho,K)
        x1=np.abs(np.trace((np.kron(M1,N1+N2)+np.kron(M2,N1-N2))@state))
        
        if x1>xtop:
            xtop=x1
            #meas=[M1,M2,N1,N2]
            fstate=state
        if track==True: print(xtop)
    print(xtop)
    return xtop,fstate#,meas

def TrAAp(m):
    rho1=np.trace(m.reshape([3,2,3,2,3,2,3,2]),axis1=0,axis2=4)
    return np.trace(rho1,axis1=1,axis2=4).reshape((4,4))

def TrBBp(m):
    rho1=np.trace(m.reshape([3,2,3,2,3,2,3,2]),axis1=1,axis2=5)
    return np.trace(rho1,axis1=2,axis2=5).reshape((9,9))

def TrApBp(m):
    rho1=np.trace(m.reshape([3,2,3,2,3,2,3,2]),axis1=2,axis2=6)
    return np.trace(rho1,axis1=2,axis2=5).reshape((6,6))

def TrAB(m):
    rho1=np.trace(m.reshape([3,2,3,2,3,2,3,2]),axis1=0,axis2=4)
    return np.trace(rho1,axis1=0,axis2=3).reshape((6,6))

def Block(vals,vecs):
    return vecs@np.sign(np.diag(vals))@np.conj(vecs.T)

Id9,Id6,Id4,Id3,Id2=np.identity(9),np.identity(6),np.identity(4),np.identity(3),np.identity(2)

def SwapMbasis(F1,mode='F'):
    if mode=='F':
        F1[:,[2,3]]=F1[:,[3,2]]
        F1[:,[3,4]]=F1[:,[4,3]]
        F1[[2,3],:]=F1[[3,2],:]
        F1[[3,4],:]=F1[[4,3],:]
    else:
        F1[[3,4],:]=F1[[4,3],:]
        F1[[2,3],:]=F1[[3,2],:]
        F1[:,[3,4]]=F1[:,[4,3]]
        F1[:,[2,3]]=F1[:,[3,2]]
    return F1

def FindMs32(N1,N2,state):
    F1=SwapMbasis(TrBBp((N1+N2)@state),mode='F')
    ev11,evecs11=eigh(F1[0:4,0:4])
    ev12,evecs12=eigh(F1[4:6,4:6])
    ev13,evecs13=eigh(F1[6:8,6:8])
    m1=np.zeros((9,9),dtype=complex)
    m1[0:4,0:4],m1[4:6,4:6],m1[6:8,6:8]=Block(ev11,evecs11),Block(ev12,evecs12),Block(ev13,evecs13)
    m1[8,8]=np.sign(F1[8,8])
    
    F2=SwapMbasis(TrBBp((N1-N2)@state))
    ev21,evecs21=eigh(F2[0:4,0:4])
    ev22,evecs22=eigh(F2[4:6,4:6])
    ev23,evecs23=eigh(F2[6:8,6:8])
    m2=np.zeros((9,9),dtype=complex)
    m2[0:4,0:4],m2[4:6,4:6],m2[6:8,6:8]=Block(ev21,evecs21),Block(ev22,evecs22),Block(ev23,evecs23)
    m2[8,8]=np.sign(F2[8,8])
    
    return SwapMbasis(m1,mode='B'),SwapMbasis(m2,mode='B')

def FindNs32(M1,M2,state):
    F1=TrAAp((M1+M2)@state)
    F2=TrAAp((M1-M2)@state)
    ev1,evec1=eigh(F1)
    ev2,evec2=eigh(F2)
    return Block(ev1,evec1), Block(ev2,evec2)

def Findrho1(Bp,rho2):
    ev,evec = eigh(TrApBp(np.kron(Id6,rho2)@Bp)[0:4,0:4])
    rho1=np.zeros((6,6),dtype=complex)
    rho1[0:4,0:4] =np.outer(evec[:,-1],np.conj(evec[:,-1]))
    return rho1

def Findrho2(Bp,rho1):
    ev,evec = eigh(TrAB(np.kron(rho1,Id6)@Bp)[0:4,0:4])
    rho2=np.zeros((6,6),dtype=complex)
    rho2[0:4,0:4] =np.outer(evec[:,-1],np.conj(evec[:,-1]))
    return rho2

def ConstBp(B,K):
    O=np.zeros(np.shape(B),dtype=complex)
    for k in K:
        O+=np.conj(k.T)@B@k
    return O

keys=[a+b+ap+bp for a in ['0','1','2'] for b in ['0','1'] for ap in ['0','1','2'] for bp in ['0','1']]

def Pad(Key,m):
    if Key=='M':
        O=np.zeros((36,36),dtype=complex)
        for i in range(36):
            for j in range(36):
                srow,scol=keys[i],keys[j]
                if srow[1]==scol[1] and srow[3]==scol[3]:
                    O[i,j]=m[int(srow[0]+srow[2],3),int(scol[0]+scol[2],3)]
    elif Key=='N':
        O=np.zeros((36,36),dtype=complex)
        for i in range(36):
            for j in range(36):
                srow,scol=keys[i],keys[j]
                if srow[0]==scol[0] and srow[2]==scol[2]:
                    O[i,j]=m[int(srow[1]+srow[3],2),int(scol[1]+scol[3],2)]
    return O

small=np.array([np.sqrt(0.99),0,0,np.sqrt(1-0.99),0,0])
rhoSmall=np.outer(small,small)

def Jump(rho,eps,pJump,state=rhoGHZ):
    p=np.random.uniform(0,1)
    if p<pJump:
        return (eps)*rho+(1-eps)*rhoSmall,True
    else: return rho,False

def SeeSaw(errTol,iters,epochs,rho1,rho2,K,eps=1e-3,pJump=1/2,maxtimes=10,stateFixed=False,track=True):
    xtop=0
    for epoch in range(epochs):

        M1=np.kron(np.kron(QutritMeasRand(),Id2),np.kron(QutritMeasRand(),Id2))
        M2=np.kron(np.kron(QutritMeasRand(),Id2),np.kron(QutritMeasRand(),Id2))
        N1=np.kron(np.kron(Id3,QubitMeasRand()),np.kron(Id3,QubitMeasRand()))
        N2=np.kron(np.kron(Id3,QubitMeasRand()),np.kron(Id3,QubitMeasRand()))
        state=ChannelState(np.kron(rho1,rho2),K)
        
        maximums,times=[2.0],[0]
        cnt=0
        chshVals=[2.0,2.0]
        
        while max(times)<maxtimes:
            
            if abs(chshVals[-1]-chshVals[0])<errTol:
                if chshVals not in maximums: 
                    maximums.append(chshVals)
                    times.append(1)
                else: times[maximums.index(chshVals)]+=1
                rho1,didIt1=Jump(rho1,eps,pJump)
                rho2,didIt2=Jump(rho2,eps,pJump)
                if didIt1==True or didIt2==True: 
                    state=ChannelState(np.kron(rho1,rho2),K)
                
            for i in range(iters):
                m1,m2=FindMs32(N1,N2,state)
                M1,M2=Pad('M',m1),Pad('M',m2)
                if stateFixed==False:
                    rho1=Findrho1(ConstBp(M1@(N1+N2)+M2@(N1-N2),K),rho2)
                    state=ChannelState(np.kron(rho1,rho2),K)
                n1,n2=FindNs32(M1,M2,state)
                N1,N2=Pad('N',n1),Pad('N',n2)
                if stateFixed==False:
                    rho2=Findrho2(ConstBp(M1@(N1+N2)+M2@(N1-N2),K),rho1)
                    state=ChannelState(np.kron(rho1,rho2),K)
                    
            x1=np.abs(np.trace((M1@(N1+N2)+M2@(N1-N2))@ChannelState(np.kron(rho1,rho2),K)))
            chshVals[cnt%2]=x1
            cnt+=1
                    
            if max(chshVals)>xtop:
                xtop=max(chshVals)
                fstate=state
        if track==True: print(xtop)
    return xtop,fstate#

my_permF =[0,1,2,3,16,17,4,5,6,7,18,19,8,9,10,11,20,21,12,13,14,15]
my_permF+=[22,23,24,25,26,27,28,29,30,31,32,33,34,35]

my_permB =[0,1,2,3,6,7,8,9,12,13,14,15,18,19,20,21,4,5,10,11,16,17]
my_permB+=[22,23,24,25,26,27,28,29,30,31,32,33,34,35]

I = np.argsort(my_permF)
J = np.argsort(my_permB)

def SwitchRhobasis(rho,mode='F'):
    if mode=='F': 
        s=rho[:,I]
        return s[I,:]
    if mode=='B':
        s=rho[J,:]
        return s[:,J]

def Findrho(B):
    ev,evecs=eigh(SwitchRhobasis(B,mode='F')[0:16,0:16])
    rho=np.zeros((36,36),dtype=complex)
    rho[0:16,0:16]=np.outer(evecs[:,-1],np.conj(evecs[:,-1]))
    return SwitchRhobasis(rho,mode='B')

def JumpFull(rho,eps,pJump):
    p=np.random.uniform(0,1)
    if p<pJump:
        return (eps)*rho+(1-eps)*np.kron(rhoSmall,rhoSmall),True
    else: return rho,False

def SeeSawFullState(errTol,iters,epochs,rho,K,eps=1e-3,pJump=1/2,maxtimes=10,stateFixed=False,track=True):
    xtop=0
    for epoch in range(epochs):
        M1=np.kron(np.kron(QutritMeasRand(),Id2),np.kron(QutritMeasRand(),Id2))
        M2=np.kron(np.kron(QutritMeasRand(),Id2),np.kron(QutritMeasRand(),Id2))
        N1=np.kron(np.kron(Id3,QubitMeasRand()),np.kron(Id3,QubitMeasRand()))
        N2=np.kron(np.kron(Id3,QubitMeasRand()),np.kron(Id3,QubitMeasRand()))
        
        state=ChannelState(rho,K)
        
        maximums,times=[1.5],[0]#a list of maxima and their times of occurance [2.,2.1],[4,2]
        cnt=0
        chshVals=[1.5,1.5]
        
        while max(times)<maxtimes:
            
            if abs(chshVals[-1]-chshVals[0])<errTol:
                if chshVals not in maximums: 
                    maximums.append(chshVals)
                    times.append(1)
                else: times[maximums.index(chshVals)]+=1
                rho,didIt1=JumpFull(rho,eps,pJump)
                if didIt1==True: 
                    state=ChannelState(rho,K)
                
            for i in range(iters):
                m1,m2=FindMs32(N1,N2,state)
                M1,M2=Pad('M',m1),Pad('M',m2)
                if stateFixed==False:
                    rho=Findrho(ConstBp(M1@(N1+N2)+M2@(N1-N2),K))
                    state=ChannelState(rho,K)
                n1,n2=FindNs32(M1,M2,state)
                N1,N2=Pad('N',n1),Pad('N',n2)
                if stateFixed==False:
                    rho=Findrho(ConstBp(M1@(N1+N2)+M2@(N1-N2),K))
                    state=ChannelState(rho,K)
                    
            x1=np.abs(np.trace((M1@(N1+N2)+M2@(N1-N2))@state))
            chshVals[cnt%2]=x1
            #print(chshVals)
            cnt+=1
                    
            if max(chshVals)>xtop:
                xtop=max(chshVals)
                fstate=state
        if track==True: print(xtop)
    return xtop,fstate#

def UntamperedSeeSawFullState(iters,epochs,rho,K,stateFixed=False,track=True):
    xtop=0
    for epoch in range(epochs):
        M1=np.kron(np.kron(QutritMeasRand(),Id2),np.kron(QutritMeasRand(),Id2))
        M2=np.kron(np.kron(QutritMeasRand(),Id2),np.kron(QutritMeasRand(),Id2))
        N1=np.kron(np.kron(Id3,QubitMeasRand()),np.kron(Id3,QubitMeasRand()))
        N2=np.kron(np.kron(Id3,QubitMeasRand()),np.kron(Id3,QubitMeasRand()))
        
        state=ChannelState(rho,K)

        for i in range(iters):
            m1,m2=FindMs32(N1,N2,state)
            M1,M2=Pad('M',m1),Pad('M',m2)
            if stateFixed==False:
                rho=Findrho(ConstBp(M1@(N1+N2)+M2@(N1-N2),K))
                state=ChannelState(rho,K)
            n1,n2=FindNs32(M1,M2,state)
            N1,N2=Pad('N',n1),Pad('N',n2)
            if stateFixed==False:
                rho=Findrho(ConstBp(M1@(N1+N2)+M2@(N1-N2),K))
                state=ChannelState(rho,K)
                    
        x1=np.abs(np.trace((M1@(N1+N2)+M2@(N1-N2))@state))
                    
        if x1>xtop:
            xtop=x1
            fstate=rho
        if track==True: print(xtop)
    return xtop,fstate#
