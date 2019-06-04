#include <vector>

void f(std::vector<int>& v) {
  v[0]; 
  v.at(0); // Guaranteed to throw a std::out_of_range exception
}

void PrettyFormat(int i, char* buf) {
  sprintf(buf, "%4d" i);
}


/*
sprintf?  What is the purpose of sprintf?

All animals are equal, but some animals are more equal than others

sprintf - changes an integer value to a human-readable string representation



 */
int main() {

  return 0;
 
}
