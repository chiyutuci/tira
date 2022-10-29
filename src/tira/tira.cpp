#include "tira/tira.h"

void Tira::execute_jobs() {
  mkdir("./relations", 0777);
  mkdir("./result", 0777);

  config_family_();

  create_relations_();
}

void Tira::config_family_() {

  // register symbols
  // symbols["m"] = get_symbol("m");
  symbols_["k1"] = get_symbol("k1");
  symbols_["p1"] = get_symbol("p1");
  symbols_["p2"] = get_symbol("p2");
  symbols_["p3"] = get_symbol("p3");
  parser reader(symbols_);

  // family name
  name_ = "box";

  // dimension d = 4 - 2 * eps
  dimension_ = get_symbol("d");

  // external momenta
  ext_vars_.push_back(get_symbol("p1"));
  ext_vars_.push_back(get_symbol("p2"));
  ext_vars_.push_back(get_symbol("p3"));

  // internal momenta
  int_vars_.push_back(get_symbol("k1"));

  // invariants
  // inv_var=get_symbol("m");

  // invariant set to 1
  // no such symbol

  // scalar products rules
  sp_rules_.append(reader("p1") * reader("p1") == reader("0"))
      .append(reader("p1") * reader("p2") == reader("50"))
      .append(reader("p1") * reader("p3") == reader("-1/2"))
      .append(reader("p2") * reader("p2") == reader("0"))
      .append(reader("p2") * reader("p3") == reader("-99/2"))
      .append(reader("p3") * reader("p3") == reader("0"));

  // propagators
  nprops_ = 4;
  // raw propagator input
  props_raw_.push_back({"k1", "0"});
  props_raw_.push_back({"k1+p1", "0"});
  props_raw_.push_back({"k1+p1+p2", "0"});
  props_raw_.push_back({"k1-p3", "0"});
  // transform to complete propagators
  for (size_t i = 0; i < nprops_; ++i) {
    auto prop_raw = props_raw_[i];
    ex prop = reader(prop_raw.first) * reader(prop_raw.first) -
              reader(prop_raw.second);
    props_.push_back(prop.expand());
  }

  // top level sector, depend on target integrals
  top_sector_ = 15; // [1,1,1,1]
  // integrand indices, depend on target integrals
  rmax_ = 6;
  smax_ = 2;
}

void Tira::create_relations_() {
  create_ibp_();
  create_li_();
  create_sym_();
}

void Tira::create_ibp_() {
  cout << "\n========== Generate IBP relations ==========\n" << endl;
  Timer timer;

  create_ibp_detail_(int_vars_); // k_i k_j
  create_ibp_detail_(ext_vars_); // k_i p_j

  // sort these relations by their length
  sort(ibps.begin(), ibps.end(),
       [](const EquationPtr &l, const EquationPtr &r) {
         return l->integrals.size() < r->integrals.size();
       });

  // save ibp relations to a file
  ofstream ibpfile("./relations/ibp");
  for (size_t i = 0; i < ibps.size(); ++i)
    ibps[i].get()->write_file(ibpfile, name_);

  cout << "\nThere are " << ibps.size() << " IBP relations\n" << endl;
  timer.end();
}

void Tira::create_li_() {
  cout << "\n========== Generate LI relations ===========\n" << endl;
  Timer timer;

  string li_file = "./relations/LI";

  timer.end();
}

void Tira::create_sym_() {}

// Integrate[ Diff[k_i^u] * (p_j^u * integrand) ] == 0
// Integrate[ Diff[k_i^u] * (k_j^u * integrand) ] == 0
void Tira::create_ibp_detail_(const vsy &vars) {

  // generate symbolic propagator indices {a1,...,an}, n=nprops
  sy *indices = new sy[nprops_];
  generate_symbols(indices, "a", nprops_);

  // coefficient of an integral in the relation.
  ex coeff;

  // j represents p_j or k_j
  // i represents k_i
  // v represents props[v] (Leibniz rule)
  for (size_t j = 0; j < vars.size(); ++j) {
    for (size_t i = 0; i < int_vars_.size(); ++i) {
      // init an ibp relation
      EquationPtr ibp(new Equation());

      // general cases
      for (size_t v = 0; v < nprops_; ++v) {
        // compute the coefficient w.r.t. props[v]
        coeff = diff(props_[v], int_vars_[i]) * (-indices[v]) * vars[j];
        coeff = coeff.expand().subs(sp_rules_, subs_options::algebraic);
        // generate a new integral
        if (coeff != 0) {
          RIntegralPtr integral(new RIntegral(nprops_));
          integral->coeff = coeff;
          integral->indices[v] = 1; // [1,0] actually represents [a1+1,a2]
          ibp->integrals.push_back(move(integral));
        }
      }
      // special case when k_j is k_i
      if (vars[j] == int_vars_[i]) {
        RIntegralPtr integral(new RIntegral(nprops_));
        integral->coeff = dimension_;
        ibp->integrals.push_back(move(integral));
      }

      // add this relation
      ibps.push_back(move(ibp));
    }
  }

  delete[] indices;
}