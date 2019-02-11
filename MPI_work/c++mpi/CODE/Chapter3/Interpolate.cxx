#ifndef INTERPOLATE
#define INTERPOLATE

#include <vector>

namespace Genfun {
  class InterpolatingPolynomial: public AbsFunction {
  public:
    // Contructor
    InterpolatingPolynomial();
    virtual double operator() (double argument) const;
    // Add a new point
    void addPoint(double x, double y);
    
  private:
    std::vector<std::pair<double, double> > xPoints;
    mutable std::vector<double> q;
  };
}

#endif INTERPOLATE
