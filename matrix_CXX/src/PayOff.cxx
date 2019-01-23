#include <algorithm>

// Constructor with a single strike paramter
PayOffCall::PayOffCall(const double _K) { K = _K;}
// Destructor (no need to use virtual keyword in source file)
PayOffCall::~PayOffCall() {}

double PayOffCall::operator() (const double S) const {
  return std::max(S-K, 0.0);
}

// Constructor with two strike parameters, upper and lower barrier
PayOffDoubleDigital::PayOffDoubleDigital(const double _U, const double _D) {
  U = _U;
  D = _D;
}

// Over-ridden operator() method, which turns PayOffDoubleDigital into a function object
double PayOffDoubleDigital::operator() (const double S) const {
  if (S >= D && S <= U) {
    return 1.0;
  } else {
    return 0.0;
  }
}
