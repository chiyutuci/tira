#ifndef TIRA_H_
#define TIRA_H_

#include "ginac/ginac.h"
#include "tira/integral.h"
#include "tira/tools.h"
#include <algorithm>
#include <sys/stat.h>

using namespace GiNaC;
using namespace std;

class Tira {
public:
  Tira() {}
  ~Tira() {}

public:
  void execute_jobs();

private:
  // configure the integral family
  void config_family_();

  void create_relations_();
  // integrate-by-parts relations
  void create_ibp_();
  void create_ibp_detail_(const vsy &vars);
  // Lorentz invariance relations
  void create_li_();
  // symmetries
  void create_sym_();

private:
  // family name
  string name_;
  // used symbols
  symtab symbols_;

  // time-space dimension in dimension regularization
  ex dimension_;
  // external momenta and loop momenta
  vsy ext_vars_, int_vars_;
  // kinematics invariant  e.g. s, m
  sy inv_var_;
  // scalar product rules  e.g. (p1+p2)^2 = s
  lst sp_rules_;

  // number of propagators
  int nprops_;
  // propagators with sperated momenta and mass: {momenta, mass}
  vector<pair<string, string>> props_raw_;
  // propagators: momenta^2 - mass
  vector<ex> props_;

  // the maximal sector
  int top_sector_;
  // integrand indices.
  // rmax: the sum of denominator indices
  // smax: the sum of numerator indices
  int rmax_, smax_;

  // IBP relations of this family
  vector<EquationPtr> ibps;
};

#endif // TIRA_H_