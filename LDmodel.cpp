#include "LDmodel.h"
LDmodel:LDmodel(){}
LDmodel:~LDmodel(){
  if(coeffs_){
    delete [] coeffs;
  }
        }

void initLD(double *coeff, int nc,int modelflag){
  coeffs_=new double [nc];
  for (int i=0;i<nc;i++){
    coeffs_[i]=coeff[i];
  }
  modelflag_=modelflag;

}

double CalLD(double* x, int nx){
  switch(modelflag_){
    case 1:
      return CallLD_quad_(double* x, int nx);
    case 2:
      return CallLD_lin_(double* x, int nx);
    case 3:
      return CallLD_sqrt_(double* x, int nx);
    case 4:
      return CallLD_4para_(double* x, int nx);
    case 5:
      return CallLD_map_(double* x, int nx);
  }
  
}

double LDmodel:CallLD_quad_(double* x, int nx){
  double ld;
  double mu = 0;
  mu = sqrt(pow(x[0],2)+ pow(x[1],2));
  ld =  1.0-(coeffs_[0]*(1-sqrt(1-mu*mu)))-(coeffs_[1]*pow((1-sqrt(1-mu*mu)),2));
  return ld;
}


double LDmodel:CallLD_lin_(double* x, int nx){
  double ld;

  return ld;
}

double LDmodel:CallLD_sqrt_(double* x, int nx){
  double ld;

  return ld;
}

double LDmodel:CallLD_4para_(double* x, int nx){
  double ld;

  return ld;
}

double LDmodel:CallLD_map_(double* x, int nx){
  double ld;

  return ld;
}
