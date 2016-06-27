#include "simplemodel.h"

int main(int argc, char *argv[]){
	double star_params[] = {0.3,0.3,2.1,0.95}; // grid size, limb1, limb2, grav1 (unused), flattening (unused)
	double planet_params[] = {0.1,0.15,500.0,0.05}; // grid size, impact, rprs, semimajor, obliquity
	int phase_points = 100;
	double phase[phase_points];
	for(int a = 0; a < phase_points; a++){
		phase[a] = asin(((a-(phase_points/2.0))/(phase_points/2.0))/500.0);
	}
	SimpleModel model(5000,200);
	model.SetupStar(star_params,5);
	model.SetupPlanet(planet_params,5);

	double *flux_out = (double*)calloc(phase_points,sizeof(double));
	model.RelativeFlux(phase,phase_points,flux_out,phase_points);
	printf("FLUX OUTPUT ======\n");
	for(int a=0;a<phase_points;a++){
		printf("%i %.9f\n",a,flux_out[a]);
	}
	free(flux_out);
	return 0;
}
