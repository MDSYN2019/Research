#include "Core.h"

std::istream& read_hw(std::istream& in, std::vector<double>& hw);

std::istream& Core::read_common(std::istream& in) {
  in >> n >> midterm >> final;
  return in;
}

std::istream& Core::read(std::istream& in) {
  read_common(in);
  read_hw(in, homework);
  return in;
}

std::istream& Grad::read(std::istream& in) {
  read_common(in);
  in >> thesis;
  read_hw(in, homework);
  return in;
}

std::istream& Audit::read(std::istream& in) {
  in >> n;
  return in;
}

bool compare(const Core& c1, const Core& c2) {
  return c1.name() < c2.name();
}

bool compare_Core_ptrs(const Core* cp1, const Core* cp2) {
  return compare(*cp1, *cp2);
}

std::string letter_grade(double grade) {
  static const double numbers[] = {
    97, 94, 90, 87, 84, 80, 77, 74, 70, 60, 0
  };

  static const char* const letters[] = {
    "A+", "A", "A-", "B+", "B", "B-", "C+", "C", "C-", "D", "F"
  };

  static const size_t ngrades = sizeof(numbers) / sizeof(*numbers);

  for (size_t i = 0; i < ngrades; ++i)
    if (grade >= numbers[i]) return letters[i];

  return "?\?\?";
}
