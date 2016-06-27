#include "Laramodel.h"
LaraModel::LaraModel(double fratio,double phi,double Req,  double ggraveq, double groteq):fratio_(fratio),phi_(phi),Req_(Req),ggraveq_(ggraveq),groteq_(groteq){
  //quad2d = new QuadTwoD(*this);
  pi_=3.1415926535897931;
  tol_ = 1.e-5;
  Maxiter_ = 50000;
  w_ = sqrt(groteq_/ggraveq);
}

LaraModel::~LaraModel(){
  //delete quad2d;
}

double LaraModel::Calgeff(double *x, int nx){
  double T;
  double d,z;
  double r[3];
  d = Determinant_(x[0],x[1]);
  z = Calzcoord_(x[0],x[1],d);
  r[0] = x[0]/Req_; r[1] = x[1]/Req_; r[2] = z/Req_;
  CalTeff_(r,1,3,&T,1);
  return T;
}



double LaraModel::Determinant_(double x,double y){
  double da, db, dc;
  double f2 = (1-fratio_)*(1-fratio_);
  double sinphi2 = sin(phi_)*sin(phi_);
  double cosphi2 = cos(phi_)*cos(phi_);
  double d;
  da = 4.*y*y * (1-f2)*(1-f2) * sinphi2 * cosphi2;
  db = cosphi2 * f2 + sinphi2;
  dc = (y*y * sinphi2 - Req_*Req_ + x*x) * f2 + y*y * cosphi2;
  d = da - 4*db*dc;
  if(d<0){
    d=0;
  }
  return d;
}

double LaraModel::Calzcoord_(double x, double y, double d){
  double za, zb;
  double f2 = (1-fratio_)*(1-fratio_);
  za = -2*y*(1-f2)*sin(phi_)*cos(phi_)+sqrt(d);
  zb = 2*(f2*cos(phi_)*cos(phi_)+sin(phi_)*sin(phi_));
  return za/zb;
}

void LaraModel::Rotate_(double *x, double *xnew){
  xnew[0] = x[0];
  xnew[1] = x[1]*cos(phi_)+x[2]*sin(phi_);
  xnew[2] = -x[1]*sin(phi_) + x[2]*cos(phi_);
  return;
}

double LaraModel::Caltau_(double theta, double r){
  double tau;
  tau = 1./3.*w_*w_*pow(r,3.)*pow(cos(theta),3.)+cos(theta)+log(tan(theta/2.));
  return tau;
}

double LaraModel::Solveangle_(double tau,double theta){
  double x = 0., x0 = tan(theta/2.),h;
  /*newton method*/
  double fx,dfx;
  for (int i=0; i < Maxiter_; i++){
    x = x0;
    //printf("x=%f\n",x);
    fx = AngleFx_(x,tau);
    if(fabs(fx)<1.e-5){
      return x;
    }
    dfx = AngleDfx_(x);
    h = fx/dfx;
    while((x-h)<0){
      h/=2.;
    }
    x-=h;
  }
  printf("in angle,exceed maxiter,fx=%f,x=%f,theta=%f\n",fx,x,theta/pi_);
  return x;
  /*bisec method*/
  //double x1=0.01, x2 = 0.99*pi_/2.,dx;
  //double f = AngleFx_(tan(x1/2.),tau);
  //double fmid = AngleFx_(tan(x2/2.),tau);
  //if (f*fmid >=0.0) {
  //  printf("Root not bracketed,f=%f,fmid=%f\n",f,fmid);
  //  throw("Root must be bracketed fro bisection in rtbis");
  //}
  //double rtb;
  //rtb = f<0.0?(dx=x2-x1,x1):(dx=x1-x2,x2);
  //for (int j=0; j< Maxiter_; j++){
  //  double xmid = rtb+dx*0.5;
  //  dx*=0.5;
  //  fmid = AngleFx_(tan(xmid/2.),tau);
  //  if(fmid<=0.) rtb = xmid;
  //  if(fabs(dx) < tol_ || fmid ==0) return rtb;
  //}
  //throw("Too many bisections in rtbis");
}
double LaraModel::AngleFx_(double x,double tau){
  double fx;
  fx = (1.-x*x)/(1.+x*x)+log(x)-tau;   
  return fx;
}
double LaraModel::AngleDfx_(double x){
  double dfx; 
  dfx = -4.*x/pow((1.+x*x),2.)+1./x;
  return dfx;
}

double LaraModel::Solveradius_(double theta){
  double x = 0., x0 = 1.0,h;
  double fx,dfx;
  /*newton method*/
  for (int i=0; i < Maxiter_; i++){
    x = x0;
    //printf("x=%f\n",x);
    fx = RFx_(x,sin(theta));
    if(fabs(fx)<tol_){
      return x;
    }
    dfx = RDfx_(x,sin(theta));
    h = fx/dfx;
    while((x-h)<0){
      h/=2.;
    }
    x-=h;
  }
  printf("in radius,exceed maxiter,fx=%f,x=%f,theta=%f\n",fx,x,theta/pi_);
  return x;
} 


double LaraModel::RFx_(double x, double sintheta){
  double fx;
  fx = 0.5*pow(x,3.)*sintheta*sintheta - (1./w_/w_+0.5)*x + 1./w_/w_;
  //printf("insolver,%f %f %f\n",fx,x,sintheta);
  return fx;
}

double LaraModel::RDfx_(double x, double sintheta){
  double dfx;
  dfx = 1.5*pow(x,2.)*sintheta*sintheta - (1./w_/w_+0.5);
  return dfx;
}

void LaraModel::CalTeff_(double *r, int np, int nd, double *T, int ng){
  
 //double Tpole = pow((1/pow(1-fratio_,4.)+pow(w_,4.)*pow(1-fratio_,2.)-2*w_*w_/(1-fratio_)),1./8.)*pow((1.-w_*w_*pow((1-fratio_),3.)),-1./6.);
  double Tpole = pow(1/pow(1-fratio_,4.),1./8.)*pow(exp(2./3.*w_*w_*pow((1-fratio_),3)),0.25);
  for (int i = 0; i<np; i++){
    double absR,absRper,theta;
    double *rtemp,*rnew;
    if(T[i] == -1){
      T[i]=0.;
      continue;
    }
    rtemp = new double [3]; //x,y,z
    rnew = new double [3]; //x0,y0,z0
    for (int j=0; j<nd; j++){
      rtemp[j] =r[nd*i+j]; 
    }
    Rotate_(rtemp,rnew);
    //double absRtemp = sqrt(rtemp[0]*rtemp[0]+rtemp[1]*rtemp[1]+rtemp[2]*rtemp[2]);
    //printf("rnew=%f %f %f\n",rnew[0],rnew[1],rnew[2]);
    absR = sqrt(rnew[0]*rnew[0]+rnew[1]*rnew[1]+rnew[2]*rnew[2]);
    //absRper = sqrt(rnew[0]*rnew[0] + rnew[2]*rnew[2]);
    //printf("absR=%f,absRper=%f\n",absR,absRtemp);
    theta = acos(rnew[1]/absR);
    //printf("theta=%f\n",theta);
    if (theta<0){
      theta=fabs(theta);
    }
    
    if (theta>pi_/2.){
      theta = pi_-theta;
    }
    if(theta>pi_/2. || theta <0){
      //printf("theta=%f\n",theta);
    }

   double tau,Thalf, tanTheta, Ta,Tb;
  if(fabs(theta)<1.e-3){
      //T[i] = pow(exp(2./3.*w_*w_*absR*absR*absR),0.25)/Tpole;
      T[i] = 1.0;
      //printf("%f %f\n",T[i],w_);
      continue;
    }
    if(fabs(theta-pi_/2.)<1.e-2)
    {
      T[i] = sqrt(2./(2.+w_*w_))*pow((1-w_*w_),1./12.)*exp(-4./3.*w_*w_/pow((2+w_*w_),3));
      continue;
    }
  
    if(cos(theta)<1.e-2){
      //T[i] = sqrt(2./(2.+w_*w_))*pow((1-w_*w_),1./12.)*exp(-4./3.*w_*w_/pow((2+w_*w_),3));
      Ta = 1./pow((1-fratio_),4.) + pow(w_,4.)*(1-fratio_)*(1-fratio_)*sin(theta)*sin(theta);
      Tb = 2.*w_*w_*sin(theta)*sin(theta)/(1-fratio_);
      T[i] = pow((Ta-Tb),1./8.) * pow((1-w_*w_*pow((1-fratio_),3)),1./6.)/Tpole;
      continue;
    }
    //inversR = pow(cos(theta[i]),2.)/pow((1-fratio_),2)+pow(sin(theta[i]),2.);
    absR = Solveradius_(theta);

   
  
   
   tau = Caltau_(theta,absR);

   Thalf = Solveangle_(tau,theta);
   tanTheta = (2*Thalf)/(1-Thalf*Thalf);
   Ta = 1./pow(absR,4.) + pow(w_,4.)*absR*absR*sin(theta)*sin(theta);
   Tb = 2.*w_*w_*sin(theta)*sin(theta)/absR;
    T[i] = pow((Ta-Tb),1./8.) * sqrt(fabs(tanTheta/tan(theta)))/Tpole;
   //T[i] = pow((Ta-Tb),1./8.) * sqrt(fabs(tanTheta/tan(theta)));
  
   //printf("%f %f\n",theta,T[i]);
    delete [] rtemp;
    delete [] rnew;
  }
}

double LaraModel::Feff_(){
  double feff, f2=(1-fratio_)*(1-fratio_);
  feff = 1- sqrt(f2*cos(phi_)*cos(phi_)+sin(phi_)*sin(phi_));
  return feff;
}
