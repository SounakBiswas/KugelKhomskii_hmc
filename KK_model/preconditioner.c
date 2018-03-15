#include "global.h"
#include "math.h"
#include "mkl.h"
#include <stdlib.h>

void calc_preconditioner(){
  int i,block;
  int link1,link2;
  double temp1,temp2;

  for(i=0;i<nsites;i++){
    link1=(i+1)%lx;
    link2=(lx+i-1)%lx;
    temp1=cosh(0.5* dtau* x_hs[link1])*cosh(0.5* dtau* x_hs[link2]);
    D_pcg[i]=alpha_pcg+temp1*temp1;
    L_pcg[i]=temp1/(1.0* D_pcg[i]);
  }

  for(block=1;block<M-1;block++ ){
    for(i=0;i<nsites;i++){
      link1=(i+1)%lx;
      link2=(lx+i-1)%lx;
      temp1=cosh(0.5* dtau* x_hs[block*nsites+link1])*cosh(0.5* dtau* x_hs[block*nsites+link2]);
      temp2=cosh(0.5* dtau* x_hs[(block-1)*nsites+link1])*cosh(0.5* dtau* x_hs[(block-1)*nsites+link2]);
      D_pcg[block *nsites+i]= alpha_pcg+temp1*temp1-temp2*temp2/(1.0 *D_pcg[(block-1)*nsites+i]);
      L_pcg[block*nsites+i]= temp1/(1.0* D_pcg[block*nsites +i]);
    }
  }

  for(i=0;i<nsites;i++){
    link1=(i+1)%lx;
    link2=(lx+i-1)%lx;
    temp1=cosh(0.5* dtau* x_hs[(M-1)*nsites+link1])*cosh(0.5* dtau* x_hs[(M-1)*nsites+link2]);
    temp2=cosh(0.5* dtau* x_hs[(M-2)*nsites+link1])*cosh(0.5* dtau* x_hs[(M-2)*nsites+link2]);
    D_pcg[(M-1)*nsites+i]=alpha_pcg+temp1*temp1 - temp2*temp2/(1.0 *D_pcg[(M-2)*nsites+i]) - temp1*temp1/(1.0 *D_pcg[i]);
    L_pcg[(M-1)*nsites+i]=temp1/(1.0* D_pcg[i]);
  }

}
void apply_preconditioner(dcomplex *x,dcomplex *y){
  double sum;
  int i,block;
  calc_preconditioner();
  //apply LT**-1 X
  for(i=0; i<nsites; i++){
    aux_pcg[i].real=x[i].real;
    aux_pcg[i].imag=x[i].imag;
  }
  for(block=1; block<M-1; block++){
    for(i=0; i<nsites; i++){
      aux_pcg[block*nsites+i].real = L_pcg[(block-1)*nsites+i]* aux_pcg[(block-1)*nsites+i].real + x[block*nsites+i].real;
      aux_pcg[block*nsites+i].imag = L_pcg[(block-1)*nsites+i]* aux_pcg[(block-1)*nsites+i].imag + x[block*nsites+i].imag;
    }
  }
  for(i=0; i<nsites; i++){
    aux_pcg[(M-1)*nsites+i].real = L_pcg[(M-2)*nsites+i]* aux_pcg[(M-2)*nsites+i].real -L_pcg[(M-1)*nsites+i]* aux_pcg[i].real+ x[(M-1)*nsites+i].real;
    aux_pcg[(M-1)*nsites+i].imag = L_pcg[(M-2)*nsites+i]* aux_pcg[(M-2)*nsites+i].imag -L_pcg[(M-1)*nsites+i]* aux_pcg[i].imag+ x[(M-1)*nsites+i].imag;
  }

  for(i=0; i<nf; i++){
    aux_pcg[i].real= aux_pcg[i].real/(D_pcg[i]);
    aux_pcg[i].imag= aux_pcg[i].imag/(D_pcg[i]);
  }


  //apply LT**-1 X
  for(i=nsites-1; i>=0; i--){
    y[(M-1)*nsites+i].real= aux_pcg[(M-1)*nsites+i].real;
    y[(M-1)*nsites+i].imag= aux_pcg[(M-1)*nsites+i].imag;
  }
  for(block=M-2; block>0; block--){
    for(i=nsites-1; i>=0; i--){
      y[block*nsites+i].real = L_pcg[block*nsites+i]* y[(block+1)*nsites+i].real + aux_pcg[block*nsites+i].real;
      y[block*nsites+i].imag = L_pcg[block*nsites+i]* y[(block+1)*nsites+i].imag + aux_pcg[block*nsites+i].imag;
    }
  }
  for(i=nsites-1; i>=0; i--){
    y[i].real = L_pcg[i]* y[nsites+i].real -L_pcg[(M-1)*nsites+i]*y[(M-1)*nsites+i].real +aux_pcg[block*nsites+i].real;
    y[i].imag = L_pcg[i]* y[nsites+i].imag -L_pcg[(M-1)*nsites+i]*y[(M-1)*nsites+i].imag +aux_pcg[block*nsites+i].imag;
  }


}
