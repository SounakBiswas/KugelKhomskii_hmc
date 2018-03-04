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

xhs=np.loadtxt("hs.dat");
vecr=np.loadtxt("vecr.dat");
veci=np.loadtxt("veci.dat");
vec=vecr+veci*1j
ctr=0;
Ka=np.zeros([ns,ns],dtype=complex)
Kb=np.zeros([ns,ns],dtype=complex)
#xhs=np.ones(nf)
#print xhs
def makemat(block) :
  K=np.zeros([ns,ns],dtype=complex)
  mat=np.zeros([ns,ns],dtype=complex)
  add=block*ns
  for i in range(0,ns):
    if(i%2==0) :
      K[i,(i+1)%ns]=-1j*xhs[add+i]
      K[(i+1)%ns,i]=+1j*xhs[add+i]
    else :
      K[i,(i+1)%ns]=-1j*xhs[add+i]
      K[(i+1)%ns,i]=+1j*xhs[add+i]
  mat=LA.expm(dtau*K*0.5)
  return mat
def makesmat(block) :
  K=np.zeros([ns,ns],dtype=complex)
  mat=np.zeros([ns,ns],dtype=complex)
  add=block*ns
  for i in range(0,ns):
    if(i%2==0) :
      Ka[i,(i+1)%ns]=-1j*xhs[add+i]
      Ka[(i+1)%ns,i]=1j*xhs[add+i]
    else :
      Kb[i,(i+1)%ns]=-1j*xhs[add+i]
      Kb[(i+1)%ns,i]=1j*xhs[add+i]
  #print "expKa"
  #print LA.expm(dtau*0.5*Ka)
  #print "expKb"
  #print LA.expm(dtau*0.5*Kb)
  #print"sparse"
  mat=np.dot(LA.expm(dtau*0.5*Ka),LA.expm(dtau*0.5*Kb))
  return mat

row=np.zeros(nelem)
col=np.zeros(nelem)
data=np.zeros(nelem,dtype=complex)
for t in range(1,M) :
  for i in range(0,ns):
    row[ctr]=t*ns+i;
    col[ctr]=t*ns+i;
    data[ctr]=1;
    ctr+=1
  mat=makesmat(t-1);
  #print "block=",t
  #print mat
  for i in range (0,ns): 
    for j in range (0,ns) : 
      row[ctr]=t*ns+i;
      col[ctr]=(t-1)*ns+j;
      data[ctr]=-mat[i,j]
      ctr+=1
t=0      
for i in range(0,ns):
  row[ctr]=t*ns+i;
  col[ctr]=t*ns+i;
  data[ctr]=1;
  ctr+=1
mat=makesmat(M-1);
for i in range (0,ns): 
  for j in range (0,ns) : 
    row[ctr]=i;
    col[ctr]=(M-1)*ns+j;
    data[ctr]=mat[i,j]
    ctr+=1
Mat=sp.coo_matrix((data,(row,col)),shape=(nf,nf),dtype=complex)
Mat=Mat.tocsr();
MatT=Mat.transpose();
MatT=MatT.conjugate();
vec2=Mat.dot(vec)
vec3=MatT.dot(vec2)
np.savetxt("vec4r.dat",np.real(vec3))
np.savetxt("vec4i.dat",np.imag(vec3))
np.savetxt("vec2r.dat",np.real(vec2))
np.savetxt("vec2i.dat",np.imag(vec2))
#vec3=(Mat.transpose()).dot(vec)
#np.savetxt("vec3.dat",vec3)
#MTM=np.dot(Mat.transpose(),Mat)
#MTMI=sLA.inv(MTM)
#vec4=(MTMI).dot(vec)
#np.savetxt("vec4.dat",vec4)
#exhs2=np.zeros(nf);
##exhs2[4]=math.exp(-dtau*root2U*xhs[4]);
#exhs2[5]=math.exp(-dtau*root2U*xhs[5]);
#
#row=np.zeros(nelem)
#col=np.zeros(nelem)
#data=np.zeros(nelem)
#ctr=0;
#for t in range(1,M) :
#  mat2=np.dot(mat,np.diag(exhs2[ns*(t-1):ns*t]))
#  for i in range (0,ns): 
#    for j in range (0,ns) : 
#      row[ctr]=t*ns+i;
#      col[ctr]=(t-1)*ns+j;
#      data[ctr]=-mat2[i,j]
#      ctr+=1
#
#t=0      
#mat2=np.dot(mat,np.diag(exhs2[ns*(M-1):]))
#for i in range (0,ns): 
#  for j in range (0,ns) : 
#    row[ctr]=t*ns+i;
#    col[ctr]=(M-1)*ns+j;
#    data[ctr]=mat2[i,j]
#    ctr+=1
#Mat=sp.coo_matrix((data,(row,col)),shape=(nf,nf),dtype='float')
#vec5=Mat.dot(vec)
#np.savetxt("vec5.dat",vec5)
#print exhs2

