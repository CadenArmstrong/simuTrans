#include "Zeipelmodel.h"
ZeipelModel::ZeipelModel(double fratio,double phi,double Req,  double ggraveq, double groteq):fratio_(fratio),phi_(phi),Req_(Req),ggraveq_(ggraveq),groteq_(groteq){
}

ZeipelModel::~ZeipelModel(){
}



double ZeipelModel::Calgeff(double *x, int nx) {
    /*x[0] = x0; x[1] = y0*/
    double g;
    double *r;
    int len = 1;
    double d,z;
    double absR,absRper,gi,gj,gz,ggrave,grote;
    double *rtemp,*rnew;
    
    r = new double [3]; /* x,y,z in the current coordinates*/
    rnew = new double [3]; /*x0,y0,z0, pole in the plane coordinates*/

    d = Determinant_(x[0],x[1]);
    if(d==-1){
      d=0;
    }
    z = Calzcoord_(x[0],x[1],d);
    r[0] = x[0]; r[1] = x[1]; r[2] = z;
    
    Rotate_(r,rnew);
    absR = sqrt(rnew[0]*rnew[0]+rnew[1]*rnew[1]+rnew[2]*rnew[2]);
    absRper = sqrt(rnew[0]*rnew[0] + rnew[2]*rnew[2]);
    ggrave = -ggraveq_*(Req_/absR)*(Req_/absR);
    grote = groteq_/(Req_/absRper);
    gi = ggrave * rnew[0]/absR + grote*rnew[0]/absRper;
    gj = ggrave * rnew[1]/absR;
    gz = ggrave * rnew[2]/absR + grote*rnew[2]/absRper;
    g = sqrt(gi*gi+gj*gj+gz*gz);
    //printf("d=%f,z=%f,gi=%f,gj=%f,gz=%f\n",d,z,gi,gj,gz);
    delete [] r;
    delete [] rnew;
    return g;
}


double ZeipelModel::Determinant_(double x,double y){
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
    d=-1;
  }
  return d;
}

inline double ZeipelModel::Calzcoord_(double x, double y, double d){
  double za, zb;
  double f2 = (1-fratio_)*(1-fratio_);
  za = -2*y*(1-f2)*sin(phi_)*cos(phi_)+sqrt(d);
  zb = 2*(f2*cos(phi_)*cos(phi_)+sin(phi_)*sin(phi_));
  return za/zb;
}

inline void ZeipelModel::Rotate_(double *x, double *xnew){
  xnew[0] = x[0];
  xnew[1] = x[1]*cos(phi_)+x[2]*sin(phi_);
  xnew[2] = -x[1]*sin(phi_) + x[2]*cos(phi_);
  return;
}



