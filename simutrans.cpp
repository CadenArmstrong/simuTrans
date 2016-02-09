#include "transitmodel.h"
#include "simplemodel.cpp"

int main(int argc, char *argv[]){
	double star_params[] = {5000,0.5,0.2,0.1,0.95}; // grid size, limb1, limb2, grav1 (unused), flattening (unused)
	double planet_params[] = {200,0.1,0.15,500.0,0.05}; // grid size, impact, rprs, semimajor, obliquity
	double phase[50];
	for(int a = 0; a < 50; a++){
		phase[a] = asin(((a-25.0)/25.0)*0.75/500.0);
	}
	SimpleModel model;
	model.SetupStar(star_params,5);
	model.SetupPlanet(planet_params,5);

	double *flux_out = (double*)calloc(50,sizeof(double));
	model.RelativeFlux(phase,50,flux_out,50);
	printf("FLUX OUTPUT ======\n");
	for(int a=0;a<50;a++){
		printf("%i %.9f\n",a,flux_out[a]);
	}
	free(flux_out);
	return 0;
}
