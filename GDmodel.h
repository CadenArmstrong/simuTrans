#ifndef GDMODEL_H_
#define GDMODEL_H_
class GDmodel {
  public:
    virtual ~GDmodel(){}
    virtual double Calgeff(double* x,int nx){};
};


#endif //GDMODEL_H_
