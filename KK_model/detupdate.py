import scipy.linalg as LA
import scipy.sparse.linalg as sLA 
import scipy.sparse as sp
import numpy as np
import math
import sys
#np.set_printoptions(precision=3)
ns=4
beta=1.0
dtau=0.1
U=1
root2U= (2*U)**0.5
M=(int)(beta/dtau)
nf=ns*M
nelem=ns*ns*M+nf

#xhs=np.random.random(nf)
#np.savetxt("xhs.dat",xhs)
xhs=np.loadtxt("xhs.dat");
ctr=0;
Ka=np.zeros([ns,ns],dtype=complex)
Kb=np.zeros([ns,ns],dtype=complex)
#xhs=np.ones(nf)
#print xhs
def makemat(block,yhs) :
  K=np.zeros([ns,ns],dtype=complex)
  mat=np.zeros([ns,ns],dtype=complex)
  add=block*ns
  for i in range(0,ns):
    if(i%2==0) :
      K[i,(i+1)%ns]=-1j*yhs[add+i]
      K[(i+1)%ns,i]=+1j*yhs[add+i]
    else :
      K[i,(i+1)%ns]=-1j*yhs[add+i]
      K[(i+1)%ns,i]=+1j*yhs[add+i]
  mat=LA.expm(dtau*K*0.5)
  return mat
def makesmat(block,yhs) :
  K=np.zeros([ns,ns],dtype=complex)
  mat=np.zeros([ns,ns],dtype=complex)
  add=block*ns
  for i in range(0,ns):
    if(i%2==0) :
      Ka[i,(i+1)%ns]=-1j*yhs[add+i]
      Ka[(i+1)%ns,i]= 1j*yhs[add+i]
    else :
      Kb[i,(i+1)%ns]=-1j*yhs[add+i]
      Kb[(i+1)%ns,i]= 1j*yhs[add+i]
  mat=np.dot(LA.expm(dtau*0.5*Kb),LA.expm(dtau*0.5*Ka))
  return mat

s2b=np.zeros([ns,2],dtype=int)
for i in range (0,ns):
  s2b[i][0]=i
  s2b[i][1]=(ns+i-2)%ns
def make_pr() :
  pr=np.zeros(ns)
  for i in range (0,ns):
    pr[i]=lam[s2b[i][0]]*lam[s2b[i][1]]
  pr[0]=pr[0]*lamp;
  pr[1]=pr[1]*lamp;
  return np.diag(pr)


A=np.eye(ns,dtype=float)
lam=np.ones(ns,dtype=float)
#lam[0]=-1;
#lam[1]=-1;
#lam[2]=-1;
lamp=1
for i in range (0,M) :
  t=makemat(i,xhs)
  A=np.dot(t,A)
pr=make_pr();
A=np.dot(pr,A)
A=A+np.eye(ns)
print LA.inv(A)
print LA.det(A)
#e,v=LA.eigh(dMDM)
#print e
#zcg(MDM,np.zeros(nf,dtype=complex),vec)
#vec2=Mat.dot(vec)
#vec3=MatD.dot(vec2)
#np.savetxt("vec4i.dat",np.imag(vec3))
#np.savetxt("vec2r.dat",np.real(vec2))
#np.savetxt("vec2i.dat",np.imag(vec2))




