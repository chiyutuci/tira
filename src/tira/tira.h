#ifndef TIRA_H_
#define TIRA_H_

#include "ginac/ginac.h"
#include "tira/tools.h"
#include <algorithm>

using namespace GiNaC;
using namespace std;

typedef possymbol sy;
typedef vector<sy> vsy;

class Tira {
public:
  Tira() {}
  ~Tira() {}

public:
  void execute_jobs();

  // configure the integral family
  void config_family();

  void create_relations();
  // integral-by-parts relations
  void create_ibp();
  // Lorentz invariance relations
  void create_li();
  // symmetries
  void create_sym();

private:
  // family name
  string name;
  // used symbols
  symtab symbols;

  // time-space dimension in dimension regularization
  ex dimension;
  // external momenta and loop momenta
  vsy ext_vars, int_vars;
  // kinematics invariant  e.g. s, m
  sy inv_var;
  // scalar product rules  e.g. (p1+p2)^2 = s
  lst sp_rules;

  // the maximal sector
  int top_sector;
  // integrand indices.
  // rmax: the sum of denominator indices
  // smax: the sum of numerator indices
  int rmax, smax;
};

#endif // TIRA_H_