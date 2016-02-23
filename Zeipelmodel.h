#ifndef ZEIPELMODEL_H
#define ZEIPELMODEL_H
#include <stdio.h>
#include <math.h>
class ZeipelModel{
  public:
    /* fratio is 1-Rpole/Req; 
     * phi is the angle of star spin to sky plane;
     * ggraveq is the gravity at equator, -GMstar/Req/Req;
     * groteq is the effective gravity at equator due to 
     * rotation at equator, Omega*Omega*Req;
     * ggraveq and grateq can be in any units, as long as 
     * they are consistent with each other, and the units Req 
     * is in.
     * */
    ZeipelModel(double fratio,double phi,double Req,  double ggraveq, double groteq);
    ~ZeipelModel();
    double Calgeff(double *x, int nx);
  private:
    double Determinant_(double x,double y);
    double Calzcoord_(double x,double y,double d);
    void Rotate_(double *x, double *xnew);
    double fratio_, phi_, Req_, ggraveq_, groteq_;
};
#endif
