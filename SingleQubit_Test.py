# -*- coding: utf-8 -*-
"""
@author: Rodrigo Araiza
"""

import Superactivation_Utils as SU
import numpy as np

'''These are tests for the function SeeSawPair'''

'''First test, setting the lost channel to p=0, we should recover 2sqrt(2)'''
#K=SU.Lost_Kraus(lam=0,where=1)
#x=SU.SeeSawPair(20,20,K,track=True)
'''Test passed'''

'''As we increase the lost channel, we should see monotonically decreasing 
chsh values. We test that at 1-1/sqrt(2) there is still violation'''
#ps=np.linspace(0,1-1/np.sqrt(2),11)
#for p in ps:
#    K=SU.Lost_Kraus(lam=p,where=1)
#    x=SU.SeeSawPair(20,50,K,track=False)
'''Test passed'''

'''Numerical seach for the CHSH bound on the loss channel'''
chshVals=[]
ps=np.linspace(0.3413073900513789,0.34130739005137906,21)
for p in ps:
    K=SU.Loss_Kraus(lam=p,where=1)
    x=SU.SeeSawPair(20,40,K,track=False,stateFixed=False)
    if x<2.0: break
    chshVals.append(x)

'''For p>=0.34130739005137906 we can't get any CHSH violation '''
