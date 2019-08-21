/* Write a function that looks like it takes anotbher functon
   as a parameter, the compiler quietly translates the parameter
   to be a pointer to a function instead.

   

 */
double (*analysis) (const std::vector<Student_info>&);
typedef double (*analysis_fp) (const std::vector<Student_info>&);

// Then we can use that type to declare our function
analysis_fp get_analysis_ptr();


bool is_negative(int n) { return n < 0;}


std::vector<int>::iterator i = find_if(v.begin(), v.end(), is_negative);



// Function pointer work

int next(int n) {
  return n + 1;
}


int (*fp) (int);

// Get the function pointer fp 
fp = &next; 
fp = next; 

// The two above statements are equialent

template<class In, class Pred>
In find_if(In begin, In end, Pred f) {
  while (begin != end && if (*begin)) {
      ++begin;
    }
  return begin;
}

