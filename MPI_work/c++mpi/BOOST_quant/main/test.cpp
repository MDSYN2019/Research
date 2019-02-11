#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

// Map library

#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <istream>


using std::map;  using std::cout;
using std::cin;  using std::string;
using std::endl; using std::vector;

// read_hw istream
istream& read_hw(istream& in, vector<double>& hw) {
  if (in) {
    hw.clear(); // Make sure the vector is clean 
    // read homework grades 
    double x;
    while (in >> x) {
      hw.push_back(x);
    }
    // clear the stream so that the input will work for the next student
    in.clear();
  }
  return in;
}


//// 1 -- Student info structs 
struct Student_info {
  string name;
  double midterm, final;
  vector<double> homework;
  
  istream& read(istream&);
  double grade() const; 
};

/* 2 -- Vec class */
vector<Student_info> vs; //empty vector
vector<double> v(100); // vector with 100 elements

// obtain the names of the types used by the vector
vector<Student_info>::const_iterator b, e;
vector<Student_info>::size_type i = 0;

for (i = 0; i != vs.size(); ++i) {
  cout << vs[i].name() << endl;
 }
b = vs.begin();
e = vs.end();

// Implementing the Vec class

template <class T>
class Vec {
public:
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef size_t size_type;
  typedef T value_type;
  typedef std::ptrdiff_t difference_type;
  typedef T& reference;
  typedef const T& const_reference;
 
  // interface
  Vec() {create();}
  explicit Vec(std::size_t n, const T& val = T()) {
    create(n,val);
  }
  // remaining interface  
private:
  // implementation
  iterator data; // first element 
  iterator limit; // last element 
};



// Scoping the members of the Student_info struct
istream& Student_info::read(istream& in) {
  in >> name >> midterm >> final; // read the name, midterm and final grades
  read_hw(in,homework);
  return in;
}

double Student_info::grade() const {
  return ::grade(midterm, final, homework);
}


vector<string> split(const string& s) {
  vector<string> ret;
  typedef string::size_type string_size;
  string_size i = 0;
  while (i != s.size())  {
    // ignore leading blanks
    // invarant: characters in range [originial i, urrent i are all spaces]
    while (i != s.size() && isspace(s[i])) 
      ++i;
    // find the end of next word
    string_size j = i;
    // invariant: none of the charactee in range
    while (j != s.size() && !isspace(s[j]))
      ++j;
    if (i != j)  {
      // copy from  starting at i and taking j-i chars
      ret.push_back(s.substr(i,j-1));
      i = j;   
    }
  }  
  return ret; 
}


map<string, vector<int> > xref(instream& in, vector<string> find_words(const strings&) = split) {
  
  string line;
  int line_number = 0;
  map<string, vector<int> > ret;
  
  while (getline(in, line)) {
    ++line_number;
    
    // break the input the words
    vector<string> words = findwords
      }  
}


int main (void) {

  
  double x = 5.0;
  double y = gsl_sf_bessel_J0(x);
  cout << "J0(%g) = %.18e" <<  x <<  y;
  
  // Counting words, using associative containers   
  string s;
  map<string, int> counters; // store each word and an associated counter
  while (cin >> s) {
    ++counters[s];
  }
  for (map<string, int>::const_iterator it = counters.begin(); it != counters.end(); ++it) {
    cout << it->first << "\t" << it->second << endl;
  }
   
  // Find all the lines that refer to each word in the input
  
  return 0;
}
