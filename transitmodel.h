#ifndef TRANSITMODEL_H
#define TRANSITMODEL_H
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_sf_ellint.h>
#define pi 3.1415926535897932384626433832795028841971693993751

// SetupStar input parameter keys
// **********************
// Model 1
// **********************
#define KEY_SS_GRID_SIZE 0 // Star setup -- number of grid points per row
#define KEY_SS_LIMB_DARKENING_1 1 // Star setup -- stellar limb darkening 1
#define KEY_SS_LIMB_DARKENING_2 2 // Star setup -- stellar limb darkening 2
#define KEY_SS_GRAVITY_DARKENING_1 3 // Star setup -- stellar gravity darkening profile
#define KEY_SS_FLATTENING 4 // Star setup -- stellar flattening
// **********************


// PlanetSetup input parameter keys
// **********************
// **********************
#define KEY_PS_GRID_SIZE 0 // Planet setup -- number of points in grid (TODO?)
#define KEY_PS_IMPACT 1 // Planet setup -- impact parameter
#define KEY_PS_RPRS 2 // Planet setup -- Radius of planet / radius of star
#define KEY_PS_SEMI_MAJOR_AXIS 3 // Planet setup -- semi major axis in R star
#define KEY_PS_OBLIQUITY 4 // Planet setup -- obliquity
#define KEY_PS_ECCENTRICITY 5 // Planet setup -- eccentricity
// **********************


using namespace std;
class TransitModel{
  public:
    TransitModel(){};
    ~TransitModel(){};
    //stellar shape, limb darkening profile, gravity darkening profile
    void SetupStar(double *star_params, int np);
    //planet shape, orbit period, impact parameter/inclination,eccentricity
    void SetupPlanet(double *planet_params, int np);
    void RelativeFlux(double *phase, int np, double *deficit_flux, int nf);
  private:
    double *star_flux_map; // 2-d flux map
    double *planet_oppacity_map; // 2-d oppacity map
    int planet_grid_size;
    int star_grid_size;
    double planet_radius; // in Rp/Rs
};
#endif
