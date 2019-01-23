#ifndef __PAYOFF__
#define __PAYOFF__

class PayOff {
public:
  PayOff(); // Default (no parameter) constructor 
  virtual ~PayOff() {};  // Virtual destructor 
  virtual double operator() (const double S) const = 0; // Pure virtual method - operator() allows the class toe called like a function 
};

class PayOffCall : public PayOff {
public:
  PayOffCall(const double K_) {};
  virtual ~PayOffCall() {};
private:
  double K;
};

class PayOffDoubleDigital : public PayOff {
private:
  double U;
  double D;
public:
  PayOffDoubleDigital(const double U_, const double D_);
  virtual ~PayOffDoubleDigital();
  // Pay-off is 1 if spot within strike barriers, 0 otherwise
  virtual double operator() (const double S) const; 
};

#endif 

  
