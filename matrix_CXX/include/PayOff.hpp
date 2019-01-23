#ifndef __PAYOFF__
#define __PAYOFF__

/*

Inheritance is a powerful concept in object-orientated programming (OOP).
It allow us to model a type of relationship between object known as is-a.

*/

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

  /*
    Virtual destructors

    All of our PayOff class and subclass and subclass destrutors have so 
    far been set to virtual, using the prefixed virtual keyword.
    
    In simple terms, a virtual destructor ensures that when derived 
    subclasses go our of scope, 
    
    
   */

  virtual ~PayOffDoubleDigital();
  // Pay-off is 1 if spot within strike barriers, 0 otherwise
  virtual double operator() (const double S) const; 
};

#endif 

  
