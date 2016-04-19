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
			this->star_flattening = 1.0-sqrt(pow(1.0-star_params[KEY_SS_FLATTENING],2)*pow(cos(star_params[KEY_SS_OBLIQUITY]),2) + pow(sin(star_params[KEY_SS_OBLIQUITY]),2)); // Effective flattening
			this->star_obliquity = star_params[KEY_SS_OBLIQUITY];
			this->max_brightness = 0.0; // IMG
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
            if (star_params[KEY_SS_G_DARK]>0){
						vec[0] = 1.0*(x-this->star_grid_size_half)/this->star_grid_size_half;
						vec[1] = 1.0*(y-this->star_grid_size_half)/this->star_grid_size_half;
            //printf("vec:%f %f\n",vec[0],vec[1]);
						GEFF = zeipel.Calgeff(vec,2);
						BB = pow(GEFF,4.*star_params[KEY_SS_G_DARK]);
            } else{
              BB=1.0; 
            }
						LD =  1.0-(star_params[KEY_SS_LIMB_DARKENING_1]*(1-sqrt(1-mu*mu)))-(star_params[KEY_SS_LIMB_DARKENING_2]*pow((1-sqrt(1-mu*mu)),2));
						if (LD*BB > this->max_brightness){ //IMG
							this->max_brightness = LD*BB; //IMG
						} // IMG
						this->star_flux_map[x+y*this->star_grid_size] = LD*BB;
						//this->star_flux_map[x+y*this->star_grid_size] = BB;
						//this->star_flux_map[x+y*this->star_grid_size] = LD;
						total_flux += star_flux_map[x+y*this->star_grid_size]*this->star_pixel_size;
					}else{
						this->star_flux_map[x+y*this->star_grid_size] =  0.;
					}
					//printf("%f ",this->star_flux_map[x+y*this->star_grid_size]);
					//printf("%d %d %f\n",x,y,this->star_flux_map[x+y*this->star_grid_size]);
					//printf("%f ",BB);
				}
				//printf("\n");
			}
			this->star_total_flux = total_flux;
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
					theta = atan2(y-this->planet_grid_size_half, x-this->planet_grid_size_half);
					mu = sqrt(pow(x-this->planet_grid_size_half,2) + pow(y-this->planet_grid_size_half,2))/this->planet_grid_size_half;
					mu = mu / (r_b/sqrt(pow(r_b*cos(theta),2)+pow(sin(theta),2)));
					if(mu <= 1.0){
						this->planet_oppacity_map[x+y*this->planet_grid_size] =  1.0;
					}else{
						this->planet_oppacity_map[x+y*this->planet_grid_size] =  0.0;
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
        if (fabs(phase[a])>pi/2. && fabs(phase[a])<1.5*pi){
          flux_out[a] = 1.0;
          continue;
        }

				planet_position_x = this->semi_major*sin(phase[a])*cos(this->obliquity); // Centre of planet x
				planet_position_y = ((1.-this->star_flattening)*this->impact_parameter) + (sin(this->obliquity)* this->semi_major*sin(phase[a])); // Centre of planet y
				theta = atan2(planet_position_y, planet_position_x); // Angle between planet centre and star centre
        //printf("%d %f %f %f\n",a,planet_position_x,planet_position_y,theta);
				float d_to_p = sqrt(pow(planet_position_x,2)+pow(planet_position_y,2)); // Distance to planet centre from star centre
				float r_of_p = (this->rp_rs)*(r_b_p)/sqrt(pow(this->rp_rs*r_b_p*cos(theta),2)+pow(this->rp_rs*sin(theta),2));// a*b/sqrt((bcost)^2 + (asint)^2)
				float r_of_s = (1.0)*(r_b_s)/sqrt(pow(r_b_s*cos(theta),2)+pow(1.0*sin(theta),2));// a*b/sqrt((bcost)^2 + (asint)^2)
        

        //printf("%d %f %f %f %f %f\n",a,planet_position_x,planet_position_y,d_to_p,r_of_p,r_of_s);
			  current_flux = 0;
				if(d_to_p <= r_of_p + r_of_s){
					if(a%15==0){
					//if(a>780 && a<782){
						// IMG
						FILE *image; // IMG
						char imgname[100];// IMG
						char imgformat[] = "trans_%.5i.ppm";// IMG
						sprintf(imgname,imgformat,a);// IMG
						image = fopen(imgname, "w+");// IMG
						int linecount = 0;// IMG
						fprintf(image, "P3\n");// IMG
						fprintf(image, "# %s\n",imgname);// IMG
						fprintf(image, "%i %i\n", this->star_grid_size, this->star_grid_size);// IMG
						fprintf(image, "100\n");// IMG
						int img_bit;// IMG
						int i,j;
						float star_space_x = 0;
						float star_space_y = 0;
						int planet_grid_x = 0;
						int planet_grid_y = 0;

						for(i = 0;i<this->star_grid_size;i++){
							star_space_x = (1.0*(i-this->star_grid_size_half))/(1.0*this->star_grid_size_half);
							for(j=0;j<this->star_grid_size;j++){
								star_space_y = (1.0*(j-this->star_grid_size_half))/(1.0*this->star_grid_size_half);
								float planet_space_x = (star_space_x - planet_position_x)/this->rp_rs;
								float planet_space_y = (star_space_y - planet_position_y)/this->rp_rs;
								planet_grid_x = (int)(this->planet_grid_size_half*(1.0+planet_space_x));
								planet_grid_y = (int)(this->planet_grid_size_half*(1.0+planet_space_y));
								img_bit = (int)(100*(this->star_flux_map[i + this->star_grid_size*j]/this->max_brightness));
								if(planet_grid_x >= 0 && planet_grid_x < this->planet_grid_size){
									if(planet_grid_y >= 0 && planet_grid_y < this->planet_grid_size){
										img_bit = img_bit*(1.0-this->planet_oppacity_map[planet_grid_x + this->planet_grid_size*planet_grid_y]);
									}

								}
								fprintf(image, " %i %i %i ", 0, img_bit,img_bit);
								linecount +=1;
								if(linecount == 4){
									fprintf(image, "\n");
									linecount = 0;
								}
							}

						}
						fclose(image);
          }
					}



				if(d_to_p <= r_of_p + r_of_s){
        
					for(int x=0;x<this->planet_grid_size;x++){
						for(int y=0;y<this->planet_grid_size;y++){
              	//star_x = ((double(x-this->planet_grid_size_half)/this->planet_grid_size_half)*this->rp_rs+planet_position_x)*this->star_grid_size_half+this->star_grid_size_half;
                //star_y = ((double(y-this->planet_grid_size_half)/this->planet_grid_size_half)*this->rp_rs+planet_position_y)*this->star_grid_size_half+this->star_grid_size_half;

							float planet_space_x = (1.0*(x-this->planet_grid_size_half)/(this->planet_grid_size_half));
							float planet_space_y = (1.0*(y-this->planet_grid_size_half)/(this->planet_grid_size_half));
							float star_space_x = (planet_space_x*this->rp_rs + planet_position_x);
							float star_space_y = (planet_space_y*this->rp_rs + planet_position_y);
							star_x = (int)round(this->star_grid_size_half*(1.0+star_space_x));
							star_y = (int)round(this->star_grid_size_half*(1.0+star_space_y));
              //if(a>897 && a<898){
              //printf("%d %d %d %d %f %f %f %f %d %d\n",x,y,star_x,star_y,planet_space_x,planet_space_y,star_space_x,star_space_y,star_x0,star_y0);
        	    //}
        //printf("%d %f %f %f\n",a,planet_position_x,planet_position_y,theta);
				//if(sqrt(pow(planet_position_x,2)+pow(planet_position_y,2)) < (1.0*r_b_s/sqrt(pow(r_b_s*cos(theta),2)+pow(1.0*sin(theta),2))) + (this->rp_rs*(this->rp_rs*r_b_p)/sqrt(pow(this->rp_rs*r_b_p*cos(theta),2)+pow(this->rp_rs*sin(theta),2)))){

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



