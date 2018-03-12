import math
import numpy as np
from numpy import linalg as LA

t=1.0;
U=1.0;
B=1
nsites=4;
lx=nsites;
ly=0;
n=2*nsites;
states=2**n;
def next_nbr(ii,xx) :
  return  (xx+ii+lx)%lx ;

def localstate(a,i) :
  if(((2**i)&a)==0) :
    return 0 
  else :
    return 1

def flip(n,i,j) :
 ret=n;
 if((n&(2**i))==0) :
   ret= ret+2**i;
 else :
   ret= ret-2**i;
 if((n&(2**j))==0) :
   ret= ret+2**j;
 else :
   ret= ret-2**j;
 return ret

def fermi_sign(a,i,j) :
  k=j;
  sign=1;
  if(i>j) :
    for k in range(j+1,i) :
      sign*=(-1)**(localstate(a,k));
  else :
    for k in range (i+1,j) :
      sign*=(-1)**(localstate(a,k));
  return sign;
 





def makeH() :
  H=np.zeros([states,states],dtype=complex)
  for i in range(0,nsites) :
    temp=np.zeros([states,states],dtype='complex')
    for a in range(0,states) :
      u=2*i;
      d=2*i+1;
      occ_u=localstate(a,u);
      occ_d=localstate(a,d);
      neigh=next_nbr(i,1);
      if localstate(a,2*i)!=localstate(a,2*neigh) :
        b=flip(a,2*i,2*neigh);
	if (localstate(a,2*i)==1) :
           temp[a,b]=1j*fermi_sign(a,2*neigh,2*i);
	else :
           temp[a,b]=-1j*fermi_sign(a,2*i,2*neigh);
      if localstate(a,2*i+1)!=localstate(a,2*neigh+1) :
        b=flip(a,2*i+1,2*neigh+1);
	if (localstate(a,2*i+1)==1) :
           temp[a,b]=1j*fermi_sign(a,2*neigh+1,2*i+1);
	else :
           temp[a,b]=-1j*fermi_sign(a,2*i+1,2*neigh+1);
    temp=np.eye(states)-0.125*np.dot(temp,temp)
    H=H+temp
  return H


def makeO1(i) :
  for a in range(0,states) :
      u=2*i;
      d=2*i+1;
      neigh=next_nbr(i,1);
      occ_u=localstate(a,u);
      occ_d=localstate(a,d);
      O1[a,a]+=occ_u-occ_d;
      if(occ_u!=occ_d) :
	b=flip(a,u,d);
	O1[a,b]+=1j;

def makeO2(i) :
  for a in range(0,states) :
      u=2*i;
      d=2*i+1;
      occ_u=localstate(a,u);
      occ_d=localstate(a,d);
      O2[a,a]+=occ_u-occ_d;
   #   if(occ_u!=occ_d) :
#	b=flip(a,u,d);
#	O2[a,b]+=1;
def makeO() :
  for a in range(0,states) :
    for i in range(0,1) :
      neigh=next_nbr(i,1);
      if localstate(a,2*i)!=localstate(a,2*neigh) :
        b=flip(a,2*i,2*neigh);
	if (localstate(a,2*i)==1) :
           O[a,b]=1j*fermi_sign(a,2*neigh,2*i);
	else :
           O[a,b]=-1j*fermi_sign(a,2*i,2*neigh);
      
O1=np.zeros([2**n,2**n])
O2=np.zeros([2**n,2**n])
O=np.zeros([2**n,2**n],dtype=complex)
H=makeH()
print LA.norm(H-(H.conjugate()).transpose())
makeO1(0)
makeO2(0)
makeO();

e,v=LA.eigh(H);
print e[0:10]
O1d=np.dot(np.conjugate(np.transpose(v)),np.dot(O1,v));
O2d=np.dot(np.conjugate(np.transpose(v)),np.dot(O2,v));
Od=np.dot(np.conjugate(np.transpose(v)),np.dot(O,v));
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
  Od_ex+=expea*Od[a,a];

print  Od_ex/(Z)



