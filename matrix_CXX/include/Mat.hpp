#ifndef __SIMPLE__MATRIX__H__
#define __SIMPLE__MATRIX__H__

/*
  Generic Programming and Template Classes
  
  In order to explain how generic programming from object-orientated 
  programming (OOP) we need to review how classes work.  Recall that clas  ses allow encapsulation of data. They provide a declaration of an interface. 

  Normal Classes are bound to the data types specified for their member data. This means hat if we want to store values or objects within that class, we must specify the type of data upfront 

  -------------------
  Generic Programming
  -------------------

  One common question that arises when discussing templates is "Why use templates over normal object orientation?"

  - Generally, more errors caught at complule time and less at run-time. Thus there isn't 
  as much need for try-catch exception blocks in your code 
*/


#include <vector>

template <typename Type = double> class SimpleMatrix {
public:
  SimpleMatrix(); // Default constructor
  // Constructor specifying rows, colummns, and a default value 
  SimpleMatrix(const int& rows, const int& cols, const Type& val);
  // Copy constructor
  SimpleMatrix(const SimpleMatrix<Type>& _rhs);
  // Assignment operator overloaded
  virtual ~SimpleMatrix(); // Destructor
  // Access to the matrix values, directly, via row and column indices
  std::vector<std::vector<Type> > get_mat() const;
  Type& value(const int& row, const int& col);

private:
  std::vector<std::vector<Type> > mat; // Use a "vector of vectors to store the values
};


#endif



