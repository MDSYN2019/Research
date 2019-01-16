#ifndef __HM1210__
#define __HM1210__

namespace HM1210 { 

  double LeJ(double sigma, double epsilon, double r,double pi,double goldden,double R) 
    
  {
    double HM; 
    
    HM = (20736*pi*goldden*epsilon*pow(sigma,12)*pow(R,3)/15625)*(5*pow(R,6) + 45*pow(R,4)*pow(r,2) + 63*pow(R,2)*pow(r,4) + 15*pow(r,6))/pow(pow(r,2) - pow(R,2),9) 
      - (62208*pi*goldden*epsilon*pow(sigma,10)*pow(R,3)/21875)*(3*pow(R,4) + 14*pow(R,2)*pow(r,2) + 7*pow(r,4))/pow(pow(r,2) - pow(R,2),7);
    return HM;
  }
  
  /* Script to calculate the interaction parameters */ 
  
  double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R)
    
  { 
    double interaction; 
    
    interaction = (746496*pi*goldden*epsilon*pow(sigma,12)*pow(R,3)/15625)*(5*pow(R,6) + 27*pow(R,4)*pow(r,2) + 27*pow(R,2)*pow(r,4) + 5*pow(r,6))*r/pow(pow(r,2) - pow(R,2),10) 
      - (124416*pi*goldden*epsilon*pow(sigma,10)*pow(R,3)/3125)*(5*pow(R,4) + 14*pow(R,2)*pow(r,2) + 5*pow(r,4))*r/pow(pow(r,2) - pow(R,2),8);
    return interaction;
  }
}
  
#endif
  
