/*

------------------------
Last Updated: 27/11/2021
------------------------

------------
Useful Links: 
-------------

- https://devblogs.microsoft.com/cppblog/using-c17-parallel-algorithms-for-better-performance/  - Parallel STL 

- http://courses.cms.caltech.edu/cs11/material/cpp/donnie/cpp-ops.html - Operator overloading 

*/

#include <iostream>
#include <iomanip>
#include <vector> // std vector libraries avaliable here 
#include <map> 
#include <algorithm>
#include <iterator>
#include <numeric>
#include <vector>
#include <string>
#include <fstream>
#include <array>
#include <cctype>
#include <chrono>
#include <random>
#include <ratio>

// Boost libraries
#include <boost/lambda/lambda.hpp>
//#include "mpi.h"

/*

What is functional programming? 
------------------------------

Most of the popular programming langauges have started introducing 
features inspired by functional language programmes. Functional programming is style of programming 
that emphasizes the evaluation of expressions, rather than the execution of commands. The expressions in these languages 
are formed by using functions to combine basic values. 

To demonstrate the difference between these styles of programming, let's start 
with a simple program implemented in the imperative style, and convert it to its 
functional equivalent. 

Imagine that you want to write a function that takes a list of files 
and calculates the number of files in each. 

 Thinking imperatively, you might analyze the steps in solving the problem as follows: 

 1. Open each file

 2. Define a counter to store the number of files

 3. Read the file one character at a time, and increase the counter every time the newline character occurs 

 General benefits of functional programming: 

 - Conciseness 

 - Correctness 

 - Efficiency 


  */

 // strings
 std::string FilePath =  "/home/sang/Desktop/GIT/Research/WIP/MPI_RDKIT_WORK/test/test.txt";

 // vectors
 std::vector<std::string> FilePathVector;

 // Simple Functions
 int count_lines(const std::string& filename) {
   std::ifstream in(filename);
   return std::count(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>(), '\n');
 }

 void change_string(std::string& input) {
   // https://en.cppreference.com/w/cpp/algorithm/transform
   std::transform(input.begin(),  input.end(), input.begin(), [](unsigned char c) -> unsigned char { return std::toupper(c);});
		 
 }


 /*
 ---------------------------------
 std::transform - what does it do?
 ---------------------------------

 std::transform applies the given function to a range and stores the results in another range, keeping the original 
 elements order and beginning at d_first 

 */

 /*
 With this solution, you're no longer concerned about exactly how the counting is implemented. You are just 
 declaring that you are counting the number of newlines that appear in the given input system.  This is the main 
 idea when writing programs in the functional style - use abstractions that let you define the intent 
 instead of specifying how to do something. 

 */

 std::vector<int> count_lines_in_files_functional(const std::vector<std::string>& files) {
   std::vector<int> results;
   for (const auto& file : files) {
     results.push_back(count_lines(file));  // We no longer care about how count_lines is implemented 
   }

   return results;
 }

 std::vector<int> count_lines_in_files(const std::vector<std::string>& files, std::string condition = "0") {
   /*
     Two nested loops are used and a few variables to keep the current state 
     of the process. Although the example is simple, it has few places where you might 
     take an error. 

     Error examples: 

     - An unininitialized variable

     - An improperly updated state 

     - wrong loop condition

     The compiler will report some of these mistakes as warnings, but the mistakes that get through 
     are usually hard to find because our brains are hard-wired to ignore them.
    */

   // We are using the imperative method here

   std::vector<int> results; // Vector to store int results 
   if (condition == "0") {  
     char c = 0;
     for (const auto& file: files) {
       int line_count = 0;
       std::ifstream in(file); // Input file stream 
       while (in.get(c)) {
 	std::cout << c;
 	if (c == '\n') { // If we encounter a newline, increase the counter
 	  line_count++;
 	}
       }
       results.push_back(line_count);    
     }

   }

   else if (condition == "1") {

     std::cout << "PLACEHOLDER";

   }

   else if (condition == "2") {

     /*
       What should be done, instead of how it should be done. 
       -----------------------------------------------------
       The std::transform - traverses the items in the files collection one by one, transforms them using the 
       count_lines function, and stores the resulting values in the results vector. 

     */ 
     std::vector<int> results(files.size());
     std::transform(files.cbegin(), files.cend(), results.begin(), count_lines); // Transofrms the files, implements count_lines function on each of the elements 
     return results; 
   }
   return results;

 }

 // Calculating averages

 double average_score(const std::vector<int>& scores, std::string condition = "0") {
   int sum = 0;

   if (condition == "0") {
     for (int score: scores) {
       sum += score;
     }
     return sum / (double)scores.size();  
   }

   /*
     The STL provides a higher-order function that can sum all the items in a collection - the 
     std::accumulate algorithm. It takes the collection and the initial value for summing, and it returns the 
     sum of the initial value and all items in the collection. 
   */

   else if (condition == "1") {
     // Link - https://en.cppreference.com/w/cpp/algorithm/accumulate
     return std::accumulate(scores.cbegin(), scores.cend(), 0) / (double)scores.size();
   }    
 }

 int f(int previous_count, char c) {
   return (c != '\n') ? previous_count : previous_count + 1; 
}

int count_lines_accumulate(const std::string& s) {
  return std::accumulate(s.cbegin(), s.cend(), 0, f);
}



/*

Operator Overloading 
--------------------

One of the nice features of C++ is that you can give special meanings to operators, when they are used with user-defined classes. 
This is called operator overloading. 


1). 

*/

class Box {
private:
  double length;
  double breadth;
  double height;
public:
  double getVolume(void) {
    return length * breadth * height;
  }
  void setLength(double len) {
    length = len; 
  }
  
  void setBreadth(double bre) {
    breadth = bre;
  }
  
  void setHeight(double hei) {
    height = hei;
  }

  Box operator+(const Box& b) { // overload the + operator for the class Box which 
    Box box;
    box.length = this->length + b.length;
    box.breadth = this->breadth + b.breadth;
    box.height = this->height + b.height;
    return box;
  }
};

int main(void) { // (void) is a recommended practise in C apparently - https://www.geeksforgeeks.org/difference-int-main-int-mainvoid/
  
  std::vector<int> NEW;

  FilePathVector.push_back(FilePath);  
  NEW = count_lines_in_files(FilePathVector); 

  for(std::vector<int>::iterator it = NEW.begin(); it != NEW.end(); ++it) {
    std::cout << *it;
  }
  
  std::string NEWSTRING;
  // Testing function for capitalizing the string in NEWSTRING

  NEWSTRING = "sdlfdsf";
  change_string(NEWSTRING);  

  std::cout << NEWSTRING << std::endl;
  std::cout << average_score(NEW);


  // Operator overloading samples

  Box Box1;
  Box Box2;
  Box Box3;
  double volume = 0.0;

  Box1.setLength(12.0);
  Box2.setLength(12.0);
  Box3.setLength(12.0);

  Box1.setBreadth(12.0);
  Box2.setBreadth(12.0);
  Box3.setBreadth(12.0);


  Box1.setHeight(12.0);
  Box2.setHeight(12.0);
  Box3.setHeight(12.0);
  
  return 0;
}
