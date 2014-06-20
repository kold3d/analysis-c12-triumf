#ifndef ENERGY_LOSS_H
#define ENERGY_LOSS_H

#include <vector>
#include <TObject.h>

class EnergyLoss : public TObject {
 public:
  EnergyLoss() : TObject() {};
  EnergyLoss(const char*, double, double);

  double CalcLoss(double,double);
  double AddBack(double,double);
  double CalcRemainder(double,double);
  double CalcRange(double,double);

 private:
  double xStepSize_;
  
  std::vector<double> energy_;
  std::vector<double> dEdx_;

 public:
  ClassDef(EnergyLoss,0)
};

#endif
