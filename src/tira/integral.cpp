#include "integral.h"

// An equation int1*coeff1+int2*coeff2==0 is represented by:
// int1*(coeff1)
// int2*(coeff2)
//
// integrals are represented by:
// name[index1,index2...indexn]
void Equation::write_file(ofstream &out, const string &name) const {
  size_t nints = integrals.size();
  size_t nprops;
  if (nints)
    nprops = integrals[0].get()->indices.size();

  for (size_t i = 0; i < nints; ++i) {
    out << name << "[";
    for (size_t j = 0; j < nprops - 1; ++j)
      out << integrals[i].get()->indices[j] << ",";
    out << integrals[i].get()->indices[nprops - 1] << "]*(";
    out << integrals[i].get()->coeff << ")" << endl;
  }
  out << endl;
}