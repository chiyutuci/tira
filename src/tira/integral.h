#ifndef INTEGRAL_H_
#define INTEGRAL_H_

#include "tira/tools.h"

class Integral;
class RIntegral;
class Equation;

typedef unique_ptr<Integral> IntegralPtr;
typedef unique_ptr<RIntegral> RIntegralPtr;
typedef unique_ptr<Equation> EquationPtr;

class Integral {
public:
  Integral() {}
  explicit Integral(unsigned n) : indices(n) {}

  ~Integral() {}

  // compare two integrals' indices
  bool same_indices(const Integral &other);

public:
  // propagator indices
  vector<int> indices;
};

// RIntegral is an integral with coefficient used in relations
class RIntegral : public Integral {
public:
  RIntegral() {}
  explicit RIntegral(unsigned n) : Integral(n) {}

  ~RIntegral() {}

public:
  // coefficient of this integral in a relation
  ex coeff;
};

class Equation {
public:
  Equation() {}
  explicit Equation(unsigned n) : integrals(n) {}

  ~Equation() {}

  // collect the duplicate integrals
  void collect_integrals();

  // write an equation to a file
  void write_file(ofstream &out, const string &name) const;

public:
  // integrals in this relation
  vector<RIntegralPtr> integrals;
};

#endif // INTEGRAL_H_