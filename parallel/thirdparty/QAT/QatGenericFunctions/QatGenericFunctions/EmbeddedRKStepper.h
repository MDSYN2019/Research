//---------------------------------------------------------------------------//
// !!                                                                     !! //
//                                                                           //
//  Copyright (C) 2016 Joe Boudreau                                          //
//                                                                           //
//  This file is part of the QAT Toolkit for computational science           //
//                                                                           //
//  QAT is free software: you can redistribute it and/or modify              //
//  it under the terms of the GNU Lesser General Public License as           //
//  published by the Free Software Foundation, either version 3 of           //
//  the License, or (at your option) any later version.                      //
//                                                                           //
//  QAT is distributed in the hope that it will be useful,                   //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of           //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            //
//  GNU Lesser General Public License for more details.                      //
//                                                                           //
//  You should have received a copy of the GNU Lesser General Public         //
//  License along with QAT.  If not, see <http://www.gnu.org/licenses/>.     //
//                                                                           //
//---------------------------------------------------------------------------//

#ifndef _EmbeddedRKStepper_hh_
#define _EmbeddedRKStepper_hh_
#include "QatGenericFunctions/AdaptiveRKStepper.h"
#include "QatGenericFunctions/ExtendedButcherTableau.h"
//
// Class EmbeddedRK stepper is a class at that implements
// Fehlberg algorithms, or embedded Runge-Kutta Algorithms.  
//
namespace Genfun {

  class EmbeddedRKStepper: public AdaptiveRKStepper::EEStepper {

  public:
    
    // Constructor:
    EmbeddedRKStepper(const ExtendedButcherTableau & tableau=CashKarpXtTableau());
    
    // Destructor:
    virtual ~EmbeddedRKStepper();
    
    // Take a single step with error estimate:
    virtual void step   (const RKIntegrator::RKData       * data, 
			 const RKIntegrator::RKData::Data & sdata, 
			 RKIntegrator::RKData::Data       & ddata, 
			 std::vector<double>              & errors) const;
    // Clone self:
    virtual EmbeddedRKStepper *clone() const;

    // Return the order:
    virtual unsigned int order() const;
  private:
    ExtendedButcherTableau tableau;
  };
}
#endif
