#include "integral.h"

bool Integral::same_indices(const Integral &other) {
  size_t lthis = indices.size();
  size_t lother = other.indices.size();
  if (lthis != lother)
    return false;
  for (size_t i = 0; i < lthis; ++i)
    if (indices[i] != other.indices[i])
      return false;
  return true;
}

void Equation::collect_integrals() {
  // collect same indices
  size_t right = integrals.size();
  for (size_t i = 0; i < right; ++i)
    for (size_t j = i + 1; j < right; ++j)
      if (integrals[i]->same_indices(*integrals[j])) {
        // collect two integrals
        integrals[i]->coeff += integrals[j]->coeff;
        integrals[i]->coeff = integrals[i]->coeff.expand();
        swap(integrals[j], integrals[right - 1]);
        --right;
        --j;
      }
  // delete zero coefficients
  for (size_t i = 0; i < right; ++i) {
    if (integrals[i]->coeff == 0) {
      swap(integrals[i], integrals[right - 1]);
      --right;
      --i;
    }
  }
  integrals.resize(right);
}

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