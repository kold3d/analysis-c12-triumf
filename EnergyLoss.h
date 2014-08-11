#ifndef ENERGY_LOSS_H
#define ENERGY_LOSS_H

#include <vector>
#include <TObject.h>
#include "spline.h"

class EnergyLoss : public TObject {
 public:
  EnergyLoss() : TObject() {};
  EnergyLoss(const char*, double);

  void SetConversion(double);
  double CalcLoss(double,double);
  double AddBack(double,double);
  double CalcRemainder(double,double);
  double CalcRange(double,double);
  double CompSimpSub(tk::spline,double,double,double,int);
  double CompSimpAdd(tk::spline,double,double,double,int);

 private:
  std::vector<double> energy_;
  std::vector<double> dEdx_;
  std::vector<double> dEdx_no_norm_;

 public:
  ClassDef(EnergyLoss,0)
};

#endif
