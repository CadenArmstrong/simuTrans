#include "simplemodel.h"
#include "Zeipelmodel.cpp"

SimpleModel::SimpleModel(void){}
SimpleModel::~SimpleModel(){};
void SimpleModel::SetupStar(double *star_params, int np){
			this->star_grid_size = (int)star_params[KEY_SS_GRID_SIZE];
			if(this->star_grid_size % 2 == 0){
		//		printf("Even sized grid detected, moving to odd size\n");
				this->star_grid_size += 1;
			}
			this->star_grid_size_half = (this->star_grid_size-1)/2;
			this->star_pixel_size = 1.0/(this->star_grid_size_half*this->star_grid_size_half);
			this->star_flux_map = (double*)calloc((this->star_grid_size*this->star_grid_size),sizeof(double));
			this->star_flattening = star_params[KEY_SS_FLATTENING];
			this->star_obliquity = star_params[KEY_SS_OBLIQUITY];
			float mu = 0;
			double theta = 0;
			float r_b = 1.0-this->star_flattening;
			long double total_flux = 0;
			double LD; // Limb Darkening
			double BB; // Black body function
			double GEFF;
			double vec[2];
			double ggraveq = pow((1.0-this->star_flattening),2);

			ZeipelModel zeipel(this->star_flattening,this->star_obliquity, 1.0,ggraveq,pow(1.0-this->star_flattening,2)*star_params[KEY_SS_GROTEQ]);

			for(int x=0;x<this->star_grid_size;x++){
				for(int y = 0;y<this->star_grid_size;y++){
					mu = sqrt(pow(x-this->star_grid_size_half,2)+ pow(y-this->star_grid_size_half,2))/star_grid_size_half;
					mu = mu / (1.0*r_b/sqrt(pow(r_b*cos(theta),2)+pow(1.0*sin(theta),2)));
					theta = atan2(y-this->star_grid_size_half, x-this->star_grid_size_half);
					if(mu <= 1.0 ){
						// QUADRATIC LIMB DARKENING LAW
						vec[0] = 1.0*(x-this->star_grid_size_half)/this->star_grid_size_half;
						vec[1] = 1.0*(y-this->star_grid_size_half)/this->star_grid_size_half;
            //printf("vec:%f %f\n",vec[0],vec[1]);
						GEFF = zeipel.Calgeff(vec,2);
						BB = pow(GEFF,4*star_params[KEY_SS_G_DARK]);
						LD =  1.0-(star_params[KEY_SS_LIMB_DARKENING_1]*(1-sqrt(1-mu*mu)))-(star_params[KEY_SS_LIMB_DARKENING_2]*pow((1-sqrt(1-mu*mu)),2));
						//this->star_flux_map[x+y*this->star_grid_size] = LD*BB;
						this->star_flux_map[x+y*this->star_grid_size] = BB;
						//this->star_flux_map[x+y*this->star_grid_size] = LD;
						total_flux += star_flux_map[x+y*this->star_grid_size]*this->star_pixel_size;
					}else{
						this->star_flux_map[x+y*this->star_grid_size] =  0;
					}
					printf("%f ",this->star_flux_map[x+y*this->star_grid_size]);
					//printf("%f ",BB);
				}
				printf("\n");
			}
			this->star_total_flux = (double)total_flux;
			printf("#Star setup complete\n");
		};
void SimpleModel::SetupPlanet(double *planet_params, int np){
			this->planet_grid_size = (int)planet_params[KEY_PS_GRID_SIZE];
			if(this->planet_grid_size % 2 == 0){
			//	printf("Even sized grid detected, moving to odd size\n");
				this->planet_grid_size += 1;
			}
			this->planet_oppacity_map = (double*)calloc((this->planet_grid_size*this->planet_grid_size),sizeof(double));
			this->planet_grid_size_half = (this->planet_grid_size-1)/2;
			this->rp_rs = planet_params[KEY_PS_RPRS];
			this->planet_pixel_size = this->rp_rs*this->rp_rs/(this->planet_grid_size_half*this->planet_grid_size_half);
			this->semi_major = planet_params[KEY_PS_SEMI_MAJOR_AXIS];
			this->impact_parameter = planet_params[KEY_PS_IMPACT];
			this->obliquity = planet_params[KEY_PS_OBLIQUITY];
			this->planet_flattening = planet_params[KEY_PS_FLATTENING];
			double mu;
			double theta = 0;
			float r_b = 1.0-this->planet_flattening;
			for(int x=0;x<this->planet_grid_size;x++){
				for(int y = 0;y<this->planet_grid_size;y++){
					mu = sqrt(pow(x-this->planet_grid_size_half,2) + pow(y-this->planet_grid_size_half,2))/planet_grid_size_half;
					mu = mu / (r_b/sqrt(pow(r_b*cos(theta),2)+pow(sin(theta),2)));
					// QUADRATIC LIMB DARKENING LAW
					theta = atan2(y-this->planet_grid_size_half, x-this->planet_grid_size_half);
					if(mu <= 1.0){
						this->planet_oppacity_map[x+y*this->planet_grid_size] =  1.0;
					}else{
						this->planet_oppacity_map[x+y*this->planet_grid_size] =  0;
					}
				}
			}
			printf("#Planet setup complete\n");
		};


void SimpleModel::RelativeFlux(double *phase, int np, double *flux_out, int npo){
			printf("#Starting integration\n");
			double planet_position_x = 0;
			double planet_position_y = 0;
			long double current_flux = 0;
			int star_x = 0;
			int star_y = 0;
			float r_b_s = 1.0-this->star_flattening;
			float r_b_p = 1.0-this->planet_flattening;
			double theta = 0;
			for(int a=0; a<np;a++){
				current_flux = 0;

				planet_position_x = this->semi_major*sin(phase[a]);
				planet_position_y = ((1.-this->star_flattening)*this->impact_parameter) + tan(this->obliquity)*planet_position_x;
				theta = atan2(planet_position_y, planet_position_x);
        //printf("%d %f %f %f\n",a,planet_position_x,planet_position_y,theta);
				if(sqrt(pow(planet_position_x,2)+pow(planet_position_y,2)) < (1.0*r_b_s/sqrt(pow(r_b_s*cos(theta),2)+pow(1.0*sin(theta),2))) + (this->rp_rs*(this->rp_rs*r_b_p)/sqrt(pow(this->rp_rs*r_b_p*cos(theta),2)+pow(this->rp_rs*sin(theta),2)))){
        //if(pow(planet_position_x,2)<1){
					for(int x=0;x<this->planet_grid_size;x++){
						star_x = (((x-this->planet_grid_size_half)/this->planet_grid_size_half)*this->rp_rs+planet_position_x)*this->star_grid_size_half+this->star_grid_size_half;
						for(int y=0;y<this->planet_grid_size;y++){
							star_y = (((y-this->planet_grid_size_half)/this->planet_grid_size_half)*this->rp_rs+planet_position_y)*this->star_grid_size_half+this->star_grid_size_half;
							if(star_x >= 0 && star_x < this->star_grid_size && star_y >= 0 && star_y < this->star_grid_size){
								current_flux += this->planet_oppacity_map[x+y*this->planet_grid_size]*this->star_flux_map[star_x + this->star_grid_size*star_y]*this->planet_pixel_size;
							}
						}
					}

					flux_out[a] = (this->star_total_flux - current_flux)/this->star_total_flux;
				}else{
					flux_out[a] = 1.0;
				}
			}
		printf("#Flux integration complete\n");
};



