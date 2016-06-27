#ifndef LDMODEL_H_
#define LDMODEL_H_
class LDmodel {
  public:
    LDmodel();
    ~LDmodel();
    void initLD(double *coeff, int nc,int modelflag=1);
    double CalLD(double* x,int nx);
  private
    double CalLD_quad_(double *x, int nx); //use quad law 
    double CalLD_lin_(double *x, int nx); // use linear law
    double CalLD_sqrt_(double *x, int nx); //use sqrt law
    double CalLD_4para_(double *x, int nx); //use four para law
    double CalLD_map_(double *x,int nx); //use an external map
    int modelflag_;
    double *coeffs_;//coeffs are the coefficients or the map depend on flags
};


#endif //GDMODEL_H_
