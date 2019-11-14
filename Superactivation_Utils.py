# -*- coding: utf-8 -*-
"""
@author: Rodrigo Araiza
"""
import numpy as np
from numpy.linalg import svd
from numpy.linalg import norm
from scipy.linalg import eigh as largest_eigh
from scipy.stats import unitary_group
from numpy.linalg import inv
import numpy.linalg as LA

'''This script is a collection of functions meant to aid in the search 
of channels that provide superactivation of CHSH-nonlocality'''

def WVec(n,d=2):
    '''Produced the W state ket for n parties and local dimension d'''
    W = np.zeros(d**n)
    for i in range(n): W[d**i] = 1
    return W/np.sqrt(n)

def GHZVec(n,d=2):
    '''Produces the GHZ state ket for n parties and local dimension d'''
    GHZ = np.zeros(d**n)
    GHZ[0]=1
    cnt = 0
    for i in range(n):
        cnt+= d*i
    GHZ[cnt]=1
    return GHZ/np.sqrt(2)

def Vac(n,d=3):
    '''Produces the vacuum density matrix'''
    Vac = np.zeros((d**n,d**n))
    Vac[0,0]=1
    return Vac

'''Partial trace functions'''
def PTr(m,subs, dims,d=3):
    '''This function takes the partial trace over one subsystem callled subs
    given the local dimensions dims'''
    #print(m.toden)
    rho_tensor=m.reshape(dims+dims)
    m_new = np.trace(rho_tensor,axis1=subs[0],axis2=(subs[0]+len(dims))).reshape((d**(len(dims)-1),d**(len(dims)-1)))
    return m_new

def TrAB(m,d=2):
    '''This function returns the partial trace over the first 2 subsystems'''
    m1=m.reshape([d,d,d,d,d,d,d,d])
    m2 = np.trace(m1,axis1=0,axis2=4)
    m3 = np.trace(m2,axis1=0,axis2=3)
    return m3.reshape((d**2,d**2))

def TrApBp(m,d=2):
    '''This function returns the partial trace over the last 2 subsystems'''
    m1=m.reshape([d,d,d,d,d,d,d,d])
    m2 = np.trace(m1,axis1=3,axis2=7)
    m3 = np.trace(m2,axis1=2,axis2=5)
    return m3.reshape((d**2,d**2))

def TrC(m,d=2):
    '''Takes trace over only one subsystem, the first one'''
    m1=m.reshape([d,d,d,d,d,d])
    return np.trace(m1,axis1=0,axis2=3).reshape((d**2,d**2))


def CHSHValue(M1,M2,N1,N2,rho):
    '''Returns the CHSH value where N1 and N2 already incorporate the states'''
    return np.abs(np.trace((M1@N1+M2@N1+M1@N2-M2@N2)@rho))

'''Functions to generate swap gates for qudits'''
def toStr(n,base=2):
   convertString = "0123"
   if n < base:
      return convertString[n]
   else:
      return toStr(n//base,base) + convertString[n%base]
  
def bg(i,size,base=2):
    s=toStr(i,base=base)
    return '0'*(size-len(s))+s

def swap_bits(i,j,b):
    'swaps bits on a binary string'
    l = list(b)
    c = l[i]
    l[i]=l[j]
    l[j]=c
    return ''.join(l)

def Swap(wire1, wire2, numWires, base=2):
    'Hard coded version of swap'
    S = np.zeros((base**numWires, base**numWires))
    bins = [bg(i,size=numWires,base=base) for i in range(base**numWires)]
    for b in bins:
        S[ int(swap_bits(wire1,wire2,b),base), int(b,base)] = 1
    return S

'''Hardcoding of the sigma matrices'''
sigmax,sigmay = np.array([[0,1],[1,0]],dtype=complex), np.array([[0,-1j],[1j,0]])
sigmaz,sigma0 = np.array([[1,0],[0,-1]],dtype=complex), np.array([[1,0],[0,1]])
sigma = [sigma0,sigmax,sigmay,sigmaz]

'''Functions to produce random initial measurements'''
def RandVec(d=3):
    r = np.random.uniform(-1,1,3)
    r = r/norm(r)
    if d==3: return r[0]*sigma[5+0]+r[1]*sigma[5+1]+r[2]*sigma[5+2]
    if d==2: return r[0]*sigma[0+1]+r[1]*sigma[1+1]+r[2]*sigma[2+1]

def RandBloch(d=3):
    T,m = np.random.uniform(-1,1,(3,3)), np.zeros((d*d,d*d),dtype=complex)
    if d==3:
        for i in range(3):
            for j in range(3):
                m += T[i,j]*np.kron(sigma[5+i],sigma[5+j])
    elif d==2:
        for i in range(3):
            for j in range(3):
                m += T[i,j]*np.kron(sigma[1+i],sigma[1+j])
    return m

def RandomMeasurement(d=2):
    U=unitary_group.rvs(d**2)
    return U@np.kron(sigmaz,sigmaz)@np.conj(U.T)

'''Functions to find the states'''

def Find_rho2(B,rho1,extendible=False):
    if extendible==True: raise ValueError('No two-symmetric extendable states allowed')
    B2=TrAB(np.kron(rho1,Id4)@B)
    ev,evec = largest_eigh(B2)
    return np.outer(evec[:,-1],np.conj(evec[:,-1]))
    
def Find_rho1(B,rho2,extendible=False):
    if extendible==True: raise ValueError('No two-symmetric extendable states allowed')
    B2=TrApBp(np.kron(Id4,rho2)@B)
    ev,evec = largest_eigh(B2)
    return np.outer(evec[:,-1],np.conj(evec[:,-1]))

def IsItPositive(rho):
    return np.all(LA.eigvals(rho) > 0)

def GetVects(rho):
    '''Gets vector and correlation matrix out of a density matrix'''
    v=np.array([np.trace(np.kron(sigma[1+i],sigma[0])@rho) for i in range(3)])
    u=np.array([np.trace(np.kron(sigma[0],sigma[1+i])@rho) for i in range(3)])
    T=np.zeros((3,3),dtype=complex)
    for i in range(3):
        for j in range(3):
            T[i,j]=np.trace(np.kron(sigma[1+i],sigma[1+j])@rho)
    return v,u,T

def SudoSign(arr):
    return np.diag(np.sign(arr))

U23 = Swap(1,2,4,base=2)
Id4 = np.identity(4)

'''Functions to find the measurements'''
def FindNs(M1,M2,fstate):
    F1=TrAB(np.kron(M1+M2,Id4)@fstate)
    F2=TrAB(np.kron(M1-M2,Id4)@fstate)
    ev1, evec1 = largest_eigh(F1)
    ev2, evec2 = largest_eigh(F2)
    return evec1@SudoSign(ev1)@np.conj(evec1.T), evec2@SudoSign(ev2)@np.conj(evec2.T)

def FindMs(N1,N2,fstate):
    F1=TrApBp(np.kron(Id4,N1+N2)@fstate)
    F2=TrApBp(np.kron(Id4,N1-N2)@fstate)
    ev1, evec1 = largest_eigh(F1)
    ev2, evec2 = largest_eigh(F2)
    #print(ev1)
    return evec1@SudoSign(ev1)@np.conj(evec1.T), evec2@SudoSign(ev2)@np.conj(evec2.T)

'''Functions related to the channels'''
def Hr(Op):
    return np.conj(Op.T)

def Loss_Kraus(lam=0.5,where=1):
    '''Kraus operators for the loss channel'''
    K =[np.identity(4)*np.sqrt(1-lam)]
    if where==1:
        K+=[np.sqrt(lam)*np.kron(np.array([[1,0],[0,0]]),np.identity(2))]
        K+=[np.sqrt(lam)*np.kron(np.array([[0,1],[0,0]]),np.identity(2))]
    elif where==2:
        K+=[np.sqrt(lam)*np.kron(np.identity(2),np.array([[1,0],[0,0]]))]
        K+=[np.sqrt(lam)*np.kron(np.identity(2),np.array([[0,1],[0,0]]))]
    return K

def Depol_Kraus(lam=0.7, where=2,d=2): #0.7 is slightly under the CHSH treshold 0.706
    '''Kraus operators for the Depolarization channel'''
    if where==2:
        K =[np.kron(sigma[4*(d-2)+0],sigma[4*(d-2)+0])*(np.sqrt(1-3*lam/d**2))]
        K+=[np.kron(sigma[4*(d-2)+0],sigma[4*(d-2)+i])*np.sqrt((lam)/d**2) for i in range(1,4)]
    elif where==1:
        K =[np.kron(sigma[4*(d-2)+0],sigma[4*(d-2)+0])*(np.sqrt(1-3*lam/d**2))]
        K+=[np.kron(sigma[4*(d-2)+i],sigma[4*(d-2)+0])*np.sqrt((lam)/d**2) for i in range(1,4)]
    return K

def BitFlip_Kraus(lam=0.5,where=1,d=2):
    '''Kraus operators for the bit-flip channel'''
    if where==1:
        K = [np.sqrt(1-lam)*np.kron(sigma[4*(d-2)+0],sigma[4*(d-2)+0])]
        K+= [np.sqrt(lam)*np.kron(sigma[4*(d-2)+1],sigma[4*(d-2)+0])]
    elif where==2:
        K = [np.sqrt(1-lam)*np.kron(sigma[4*(d-2)+0],sigma[4*(d-2)+0])]
        K+= [np.sqrt(lam)*np.kron(sigma[4*(d-2)+0],sigma[4*(d-2)+1])]
    return K

def PhaseFlip_Kraus(lam=0.5,where=1,d=2):
    '''Kraus operators for the phase-flip channel'''
    if where==1:
        K = [np.sqrt(1-lam)*np.kron(sigma[4*(d-2)+0],sigma[4*(d-2)+0])]
        K+= [np.sqrt(lam)*np.kron(sigma[4*(d-2)+3],sigma[4*(d-2)+0])]
    elif where==2:
        K = [np.sqrt(1-lam)*np.kron(sigma[4*(d-2)+0],sigma[4*(d-2)+0])]
        K+= [np.sqrt(lam)*np.kron(sigma[4*(d-2)+0],sigma[4*(d-2)+3])]
    return K

def AmplitudeDamping_Kraus(lam=1/2,where=1,d=2):
    '''Kraus operators for the amplitude damping channel'''
    if where==1 and d==2:
        K = [np.kron(np.array([[1,0],[0,np.sqrt(1-lam)]]),np.identity(2))]
        K+= [np.kron(np.array([[0,np.sqrt(lam)],[0,0]]),np.identity(2))]
    elif where==2 and d==2:
        K = [np.kron(np.identity(2),np.array([[1,0],[0,np.sqrt(1-lam)]]))]
        K+= [np.kron(np.identity(2),np.array([[0,np.sqrt(lam)],[0,0]]))]
    return K

def Identity_Kraus(d=2):
    return [np.identity(d*d)]
  
def Update_State(rho1,rho2,K1,K2,base):
    '''Gets the state K.rho.K^dagger where K are Kraus operators'''
    state=np.zeros((base**4,base**4),dtype=complex)
    for k1 in range(len(K1)):
        for k2 in range(len(K2)):
                state+= np.kron(K1[k1]@rho1@Hr(K1[k1]),K2[k2]@rho2@Hr(K2[k2]))
    return state

def DoesItViolate(rho,d=2):
    '''Checks if a state rho violates the Horodecki condition for the CHSH
    value'''
    if d!=2: raise ValueError('Only qubits supported for now')
    Tm=np.zeros((3,3),dtype=complex)
    for i in range(3):
        for j in range(3):
            Tm[i,j]=np.trace(np.kron(sigma[1+i],sigma[1+j])@rho)
    U=np.conj(Tm.T)@Tm
    ev1, evec1 = largest_eigh(U, eigvals=(3-2,3-1))
    if np.abs(ev1[0]+ev1[1])<=1.0: return False
    else: return True

GHZ=np.array([1,0,0,-1])/np.sqrt(2) 
rhoGHZ=np.outer(GHZ,GHZ)

def ConstructBp(Op,KS):
    '''The matrix Bp is the effective Bell-operator due to the presence 
    of the channel'''
    O=np.zeros(np.shape(Op),dtype=complex)
    for i in range(len(KS)):
            O+= Hr(KS[i])@Op@KS[i]
    return O

def Jump(rho,eps,pJump,state=rhoGHZ):
    '''Performs step (3) in the paper's see-saw algorithm for state that is 
    of the form rhoAB'''
    p=np.random.uniform(0,1)
    if p<pJump:
        return (eps)*rho+(1-eps)*rhoGHZ,True
    else: return rho,False

def ChannelsSeeSaw(errTol,epochs,K1,K2,iters=10,eps=1e-3,pJump=1/2,maxtimes=10,s1=np.array([]),s2=np.array([]),statesFixed=False,track=True,ext=False):
    xTOP=0
    KS=[U23@np.kron(K1[i],K2[j])@U23 for i in range(len(K1)) for j in range(len(K2))]
    
    for epoch in range(epochs):
        #Initializing
        M1,M2=RandomMeasurement(),RandomMeasurement()
        N1,N2=RandomMeasurement(),RandomMeasurement()
        Bp=U23@ConstructBp(np.kron(M1,N1+N2)+np.kron(M2,N1-N2),KS)@U23
        
        if np.shape(s1)==(0,): rho1 = Find_rho1(Bp,rhoGHZ,extendible=False)
        elif np.shape(s1)!=(0,): rho1=s1
        if np.shape(s2)==(0,): rho2 = Find_rho2(Bp,rho1,extendible=False)
        elif np.shape(s2)!=(0,): rho2=s2
        state=U23@Update_State(rho1,rho2,K1,K2,2)@U23
        
        maximums,times=[2.0],[0]#a list of maxima and their times of occurance [2.,2.1],[4,2]
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
                    state=U23@Update_State(rho1,rho2,K1,K2,2)@U23
                    
            for i in range(iters):
                N1,N2=FindNs(M1,M2,state)
                Bp=U23@ConstructBp(np.kron(M1,N1+N2)+np.kron(M2,N1-N2),KS)@U23
                if statesFixed==False:
                    rho1=Find_rho1(Bp,rho2,extendible=ext)
                    state=U23@Update_State(rho1,rho2,K1,K2,2)@U23
                M1,M2=FindMs(N1,N2,state)
                Bp=U23@ConstructBp(np.kron(M1,N1+N2)+np.kron(M2,N1-N2),KS)@U23
                if statesFixed==False: 
                    rho2=Find_rho2(Bp,rho1,extendible=ext)
                    state=U23@Update_State(rho1,rho2,K1,K2,2)@U23
                Bp=U23@ConstructBp(np.kron(M1,N1+N2)+np.kron(M2,N1-N2),KS)@U23

            x1 = np.abs(np.trace((np.kron(M1,N1+N2)+np.kron(M2,N1-N2))@state))
            chshVals[cnt%2]=x1
            cnt+=1
                    
            if max(chshVals)>xTOP:
                xTOP=max(chshVals)
                states=[rho1,rho2]
        if track==True: print(xTOP)
    
    print("CHSH Value = %f" %xTOP)
    return xTOP,states

def ChannelState(psi,K1,K2):
    '''Taking the fill state, it passes it through the channel'''
    O=np.zeros(np.shape(psi),dtype=complex)
    for k1 in K1:
        for k2 in K2:
            O+= np.kron(k1,k2)@psi@np.conj(np.kron(k1,k2).T)
    return O

def FindFullState(B):
    '''Finds optimal state rho_ABA'B' based on the Bell operator's
    largest eigenvalue eigenvector'''
    ev, evec = largest_eigh(B)
    return np.outer(evec[:,-1],np.conj(evec[:,-1]))

fullGHZ=np.kron(rhoGHZ,rhoGHZ)

def FullJump(rho,eps,pJump):
    '''Performs step (3) in the paper's see-saw algorithm for state that is 
    of the form rhoABA'B'''
    p=np.random.uniform(0,1)
    if p<pJump:
        return (eps)*rho+(1-eps)*fullGHZ,True
    else: return rho,False
    
def FullStateSeeSaw(errTol,epochs,K1,K2,iters=10,eps=1e-3,pJump=1/2,maxtimes=10,track=True,ext=False):
    '''Performs the see-saw algorithm for a state of the form rhoABA'B' and not for 
    rhoAB.rhoA'B' '''
    xTOP=0
    KS=[U23@np.kron(K1[i],K2[j])@U23 for i in range(len(K1)) for j in range(len(K2))]
    
    for epoch in range(epochs):
        M1,M2=RandomMeasurement(),RandomMeasurement()
        N1,N2=RandomMeasurement(),RandomMeasurement()
        Bp=U23@ConstructBp(np.kron(M1,N1+N2)+np.kron(M2,N1-N2),KS)@U23
        rho=FindFullState(Bp)
        state=U23@ChannelState(rho,K1,K2)@U23
        
        maximums,times=[2.0],[0]
        cnt=0
        chshVals=[2.0,2.0]
        while max(times)<maxtimes:
            #print(times)
            for i in range(iters):
                N1,N2=FindNs(M1,M2,state)
                Bp=U23@ConstructBp(np.kron(M1,N1+N2)+np.kron(M2,N1-N2),KS)@U23
                rho=FindFullState(Bp)
                M1,M2=FindMs(N1,N2,state)
                state=U23@ChannelState(rho,K1,K2)@U23
                
                Bp=U23@ConstructBp(np.kron(M1,N1+N2)+np.kron(M2,N1-N2),KS)@U23
                state=U23@ChannelState(rho,K1,K2)@U23
        
            x1 = np.abs(np.trace((np.kron(M1,N1+N2)+np.kron(M2,N1-N2))@state))
            chshVals[cnt%2]=x1
            #print(chshVals)
            cnt+=1
            if abs(chshVals[-1]-chshVals[0])<errTol:
                if chshVals not in maximums: 
                    maximums.append(chshVals)
                    times.append(1)
                else: times[maximums.index(chshVals)]+=1
                rho,didIt1=FullJump(rho,eps,pJump)
                if didIt1==True: 
                    state=U23@ChannelState(rho,K1,K2)@U23
            
        x1 = np.abs(np.trace((np.kron(M1,N1+N2)+np.kron(M2,N1-N2))@state))
        if x1>xTOP: 
            xTOP=x1
            dastate=rho
        if track==True: print(xTOP)
    
    #print('See-Saw algorithm finalized')
    print("CHSH Value = %f" %xTOP)
    return xTOP,dastate

def TrABAp(m):
    '''This function returns the partial trace over the first 2 subsystems'''
    m1=m.reshape([2,2,2,2,2,2,2,2])
    m2 = np.trace(m1,axis1=0,axis2=4)
    m3 = np.trace(m2,axis1=0,axis2=3)
    m4 = np.trace(m3,axis1=0,axis2=2)
    return m4.reshape((2,2))
vac=np.array([[1,0,0,0],[0]*4,[0]*4,[0]*4])
U24=Swap(1,3,4, base=2)
U34=Swap(2,3,4, base=2)

'''Functions to be used for unital channels'''
def LamdasToKraus(ls,where=1):
    '''Builds a unital channel out of three parameters'''
    ps=1/2*np.array([[0,-1,-1],[-1,0,-1],[-1,-1,0]])@ls
    if where==1:
        K =[np.sqrt(ps[0])*np.kron(sigmax,np.identity(2))]
        K+=[np.sqrt(ps[1])*np.kron(sigmay,np.identity(2))]
        K+=[np.sqrt(ps[2])*np.kron(sigmaz,np.identity(2))]
    if where==2:
        K =[np.sqrt(ps[0])*np.kron(np.identity(2),sigmax)]
        K+=[np.sqrt(ps[1])*np.kron(np.identity(2),sigmay)]
        K+=[np.sqrt(ps[2])*np.kron(np.identity(2),sigmaz)]
    return K
        
def IsUnitalCHSH_Breaking(ls):
    '''Check if a unital channel is CHSH-breaking'''
    s=sorted(ls)
    if s[-1]**2+s[-2]**2<=1: return True
    else: return False

def Unital(rho,ls):
    v,u,T=GetVects(rho)
    s=np.identity(4,dtype=complex)
    for i in range(3):
        s+=v[i]*ls[i]*np.kron(sigma[1+i],sigma0)+u[i]*np.kron(sigma0,sigma[1+i])
    for i in range(3):
        for j in range(3):
            s+=ls[i]*T[i,j]*np.kron(sigma[1+i],sigma[1+j])
    return s/4

def Unital(M,ls1,ls2):
    v,u,T=GetVects(M)
    s=np.identity(4,dtype=complex)
    for i in range(3):
        s+=v[i]*ls1[i]*np.kron(sigma[1+i],sigma0)+u[i]*ls2[i]*np.kron(sigma0,sigma[1+i])
    for i in range(3):
        for j in range(3):
            s+=ls1[i]*ls2[j]*T[i,j]*np.kron(sigma[1+i],sigma[1+j])
    return s/4
        
def SeeSawUnital(epochs,iters,ls1,ls2,s=1,track=True,ext=False):
    '''We assume the unitaries are applied on channels A, A' '''
    xTOP=0
    base=2
    
    for epoch in range(epochs):
        #Initializing
        M1,M2=RandomMeasurement(d=base),RandomMeasurement(d=base)
        N1,N2=RandomMeasurement(d=base),RandomMeasurement(d=base)
        
        Bp=np.kron(Unitary2(M1,ls1,ls2),N1+N2)+np.kron(Unitary2(M2,ls1,ls2),N1-N2)
        
        rho1,rho2 = rhoGHZ,rhoGHZ
        
        state=U23@np.kron(Unitary(rho1,ls1),Unitary(rho2,ls2))@U23
        
        for i in range(iters):
            N1,N2=FindNs(M1,M2,state,d=base)
            M1,M2=FindMs(N1,N2,state,d=base)
            Bp=np.kron(Unitary2(M1,ls1,ls2),N1+N2)+np.kron(Unitary2(M2,ls1,ls2),N1-N2)
            
            rho1=Find_rho1(Bp,rho1,rho2,step=s,d=base,extendible=ext)
            rho2=Find_rho2(Bp,rho1,rho2,step=s,d=base,extendible=ext)
            state=U23@np.kron(Unitary(rho1,ls1),Unitary(rho2,ls2))@U23
                
        x1 = np.abs(np.trace((np.kron(M1,N1+N2)+np.kron(M2,N1-N2))@state))
        
        if x1>xTOP:
            xTOP=x1
            #states = [rho1,rho2]
            #meas   = [M1,M2,N1,N2]

        if track==True: print(xTOP)
    
    print('See-Saw algorithm finalized')
    print("CHSH Value = %f" %xTOP)
    return xTOP#,states,meas