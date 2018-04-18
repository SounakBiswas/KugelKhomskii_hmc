import math 
import numpy as np
import numpy.linalg as LA
import time
J=1.0
lx=4;
n=2*lx
states=2**n;
B=1.0
def getbit(state,i) :
    if(((2**i)&state)==0) :
        return 0 
    else :
        return 1
def getspin(state,i) :
    if(((2**i)&state)==0) :
        return -1 
    else :
        return 1
def flip(state,i,j) :
    b=2**i+2**j;
    return state^b;
state_count=0;
#Calculate Hamiltonian
H=np.zeros([states,states]);
for i in range (0,lx) :
    temp1=np.zeros([states,states])
    temp2=np.zeros([states,states])
    for state in range(0,states) :
            nbr=(i+1)%lx
            temp1[state,state]+= 1+getspin(state,2*i) * getspin(state,2*nbr) ;
            if(getbit(state,2*i)!=getbit(state,2*nbr)) :
                newstate=flip(state,2*i,2*nbr)
                temp1[state,newstate]+=1 

            temp2[state,state]+= 1.0+getspin(state,2*i+1) * getspin(state,2*nbr+1) ;
            if(getbit(state,2*i+1)!=getbit(state,2*nbr+1)) :
                newstate=flip(state,2*i+1,2*nbr+1)
                temp2[state,newstate]+=1 
    H+=np.dot(temp1,temp2)

O1=np.zeros([states,states]);
for i in range (0,1) :
    for state in range(0,states) :
            nbr=(i+2)%lx
            O1[state,state]+= getspin(state,2*i) * getspin(state,2*nbr) ;
            if(getbit(state,2*i)!=getbit(state,2*nbr)) :
                newstate=flip(state,2*i,2*nbr)
                O1[state,newstate]+=1 

#################################
print "ready"
e,v=LA.eigh(H);
Z=0
E=0
Od_ex=0
O1d=np.dot(np.conjugate(np.transpose(v)),np.dot(O1,v));
Z=0
E=0
Od_ex=0
for a in range(0,2**n) :
  Z+=math.exp(-B*e[a]);
  E+=e[a]*math.exp(-B*e[a]);
E=E/Z;
minim=0.000001/B;
for a in range(0,2**n) :
  expea=math.exp(-B*e[a]); 
  Od_ex+=expea*O1d[a,a];

print  Od_ex/(Z)

