#ifndef _VANILLA_OPTION_H
#define _VANILLA_OPTION_H

class VanillaOption {
private:
  void init();
  void copy(const VanillaOption& rhs);
  double K;
  double r;
  double T;
  double S;
  double sigma;
public:
  VanillaOption(); // Default constructor - has no parameters
  VanillaOption(const double& _K, const double& _r, const double& _T, const double& _S, const double& _sigma);
  VanillaOption(const VanillaOption& rhs); // Copy constructor 
  VanillaOption& operator=(const VanillaOption& rhs); // Assignment 
  virtual ~VanillaOption(); // Destructor 
};

#endif 
