#ifndef SIMPLEMODEL_H
#define SIMPLEMODEL_H
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
#define KEY_SS_FLATTENING 3 // Star setup -- stellar flattening
#define KEY_SS_OBLIQUITY 4 // the angle between the stellar spin with the sky plane
#define KEY_SS_GROTEQ 5 // Constant for rotational rate of star, gravitation constant, and stellar mass
#define KEY_SS_G_DARK 6 // Gravity darkening parameter
#define KEY_SS_GDFLAG 7 // Gravity darkening parameter

// **********************


// PlanetSetup input parameter keys
// **********************
// **********************
#define KEY_PS_GRID_SIZE 0 // Planet setup -- number of points in grid
#define KEY_PS_IMPACT 1 // Planet setup -- impact parameter
#define KEY_PS_RPRS 2 // Planet setup -- Radius of planet / radius of star
#define KEY_PS_SEMI_MAJOR_AXIS 3 // Planet setup -- semi major axis in R star
#define KEY_PS_OBLIQUITY 4 // Planet setup -- obliquity
#define KEY_PS_ECCENTRICITY 5 // Planet setup -- eccentricity
#define KEY_PS_FLATTENING 6
// **********************


using namespace std;
class SimpleModel{
	public:
		SimpleModel(int star_grid_size, int planet_grid_size);
		~SimpleModel();
		void SetupStar(double *star_params, int np);
		void SetupPlanet(double *planet_params, int np);
		void RelativeFlux(double *phase, int np, double *flux_out, int npo);
	private:
		double *star_flux_map;
		double *planet_oppacity_map;
		long double star_total_flux;
		double rp_rs; // Radius of planet over radius of star
		double star_pixel_size;
		double planet_pixel_size;
		double semi_major;
		double impact_parameter;
		double obliquity;
		double planet_flattening;
		double star_flattening;
		double star_obliquity;
		int star_grid_size;
		int star_grid_size_half;
		int planet_grid_size;
		int planet_grid_size_half;
		double max_brightness;

};

#endif
