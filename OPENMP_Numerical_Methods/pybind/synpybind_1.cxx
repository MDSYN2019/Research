#include <pybind11/pybind11.h>
#include <string>

#include "Eigen/Dense"
#include "Eigen/LU"
#include "Eigen/Core"

namespace py = pybind11;


//struct Pet {
  //   Pet(const std::string &name) : name(name) { }
  //    void setName(const std::string &name_) { name = name_; }
  //    const std::string &getName() const { return name; }
//    std::string name;
  //};

struct Pet {
  Pet(const std::string &name) : name(name) { }
  void setName(const std::string &name_) { name = name_; }
  const std::string &getName() const { return name; }
  std::string name;
};

struct Dog : Pet {
    Dog(const std::string &name) : Pet(name) { }
    std::string bark() const { return "woof!"; }
};

PYBIND11_MODULE(example2, m) {
    py::class_<Pet>(m, "Pet")
        .def(py::init<const std::string &>())
        .def("setName", &Pet::setName)
        .def("getName", &Pet::getName);
}


