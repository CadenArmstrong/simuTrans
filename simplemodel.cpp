#include "transitmodel.h"

class SimpleModel: public TransitModel{
	public:
		SimpleModel();
		~SimpleModel();
		void SetupStar(double *star_params, int np);
		void SetupPlanet(double *planet_params, int np);
		void RelativeFlux(double *phase, int np, double* flux_out, int npo);
	private:
		double *star_flux_map;
		double *planet_oppacity_map;
		double star_total_flux;
		double rp_rs; // Radius of planet over radius of star
		double star_pixel_size;
		double planet_pixel_size;
		double semi_major;
		double impact_parameter;
		double obliquity;
		int star_grid_size;
		int star_grid_size_half;
		int planet_grid_size;
		int planet_grid_size_half;

};

SimpleModel::SimpleModel(void){}
SimpleModel::~SimpleModel(){};
void SimpleModel::SetupStar(double *star_params, int np){
	this->star_grid_size = (int)star_params[KEY_SS_GRID_SIZE];
	this->star_pixel_size = 1.0;
	if(this->star_grid_size % 2 == 0){
		printf("Even sized grid detected, moving to odd size\n");
		this->star_grid_size += 1;
	}
	this->star_grid_size_half = (this->star_grid_size-1)/2;
	this->star_flux_map = (double*)calloc((star_params[KEY_SS_GRID_SIZE]*star_params[KEY_SS_GRID_SIZE]),sizeof(double));

	float mu = 0;
	long double total_flux = 0;
	for(int x=0;x<this->star_grid_size;x++){
		for(int y = 0;y<this->star_grid_size;y++){
			mu = sqrt(x*x + y*y)/star_grid_size;
			// QUADRATIC LIMB DARKENING LAW
			if(mu <= 1.0){
				this->star_flux_map[x+y*this->star_grid_size] =  1-star_params[KEY_SS_LIMB_DARKENING_1]*(1-mu)-star_params[KEY_SS_LIMB_DARKENING_2]*(1-mu)*(1-mu);
				total_flux += star_flux_map[x+y*this->star_grid_size]*this->star_pixel_size;
			}else{
				this->star_flux_map[x+y*this->star_grid_size] =  0;
			}
		}
	}
	this->star_total_flux = (double)total_flux;
	printf("Star setup complete\n");
};
void SimpleModel::SetupPlanet(double *planet_params, int np){
	this->planet_grid_size = (int)planet_params[KEY_PS_GRID_SIZE];
	this->planet_oppacity_map = (double*)calloc((planet_params[KEY_SS_GRID_SIZE]*planet_params[KEY_PS_GRID_SIZE]),sizeof(double));
	if(this->planet_grid_size % 2 == 0){
		printf("Even sized grid detected, moving to odd size\n");
		this->planet_grid_size += 1;
	}
	this->planet_grid_size_half = (this->planet_grid_size-1)/2;
	this->rp_rs = planet_params[KEY_PS_RPRS];
	this->planet_pixel_size = this->rp_rs*this->rp_rs;
	this->semi_major = planet_params[KEY_PS_SEMI_MAJOR_AXIS];
	this->impact_parameter = planet_params[KEY_PS_IMPACT];
	this->obliquity = planet_params[KEY_PS_OBLIQUITY];
	double mu;
	for(int x=0;x<this->planet_grid_size;x++){
		for(int y = 0;y<this->planet_grid_size;y++){
			mu = sqrt(x*x + y*y)/planet_grid_size;
			if(mu <= 1.0){
				this->planet_oppacity_map[x+y*this->planet_grid_size] =  1.0;
			}else{
				this->planet_oppacity_map[x+y*this->planet_grid_size] =  0;
			}
		}
	}
	printf("Planet setup complete\n");
};

void SimpleModel::RelativeFlux(double *phase, int np, double *flux_out, int npo){
	printf("Starting integration\n");
	double planet_position_x = 0;
	double planet_position_y = 0;
	double current_flux = 0;
	for(int a=0; a<np;a++){
		printf("Phase: %f\n",phase[a]);
		planet_position_x = this->semi_major*sin(phase[a]);
		planet_position_y = this->impact_parameter + tan(this->obliquity)*planet_position_x;
		for(int x=0;x<this->planet_grid_size;x++){
			for(int y=0;y<this->planet_grid_size;y++){
				int star_x = (((x-this->planet_grid_size_half)/this->planet_grid_size)*this->rp_rs+planet_position_x)*1.5*this->star_grid_size;
				int star_y = (((y-this->planet_grid_size_half)/this->planet_grid_size)*this->rp_rs+planet_position_y)*1.5*this->star_grid_size;
				if(star_x >= 0 && star_x < this->star_grid_size && star_y >= 0 && star_y < this->star_grid_size){
					current_flux += this->planet_pixel_size*this->planet_oppacity_map[x+y*this->planet_grid_size]*this->star_flux_map[star_x + this->star_grid_size*star_y];
				}
			}
		}
		flux_out[a] = this->star_total_flux - current_flux;
		printf("TEST\n");
	}
	printf("Flux integration complete\n");
};
