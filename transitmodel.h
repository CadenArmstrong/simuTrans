#ifndef TRANSITMODEL_H
#define TRANSITMODEL_H
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_sf_ellint.h>
#define pi 3.1415926535897932384626433832795028841971693993751
using namespace std;
class Transitmodel{
  public:
    Transitmodel();
    ~Transitmodel();
    //stellar shape, limb darkening profile, gravity darkening profile
    void Setup_star(double *star_params, int np);
    //planet shape, orbit period, impact parameter/inclination,eccentricity
    void Setup_planet(double *planet_params, int np);
    void relativeFlux(double *phase, int np, double *deficitFlux, int nf);
  private:
    
};
#endif
