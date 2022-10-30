#include "tira/tira.h"

void Tira::execute_jobs() {
  mkdir("./relations", 0777);
  mkdir("./result", 0777);

  config_family_();
  prepare_family_();

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
  // limit: must input simple sps, i.e. p1*p2==s/2, cannot input (p1+p2)^2==s
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
    props_.push_back(prop.expand().subs(sp_rules_, subs_options::algebraic));
  }
  cout << props_ << endl;

  // top level sector, depends on target integrals
  top_sector_ = 15; // [1,1,1,1]
  // integrand indices, depends on target integrals
  rmax_ = 6;
  smax_ = 2;
}

void Tira::prepare_family_() {
  size_t nint = int_vars_.size();
  size_t next = ext_vars_.size();
  size_t nint2 = nint * (nint + 1) / 2;

  // prepare: symbols of propagators dd_i
  matrix propsy(nprops_, 1);
  for (size_t i = 0; i < nprops_; ++i)
    propsy(i, 0) = get_symbol("dd" + to_string(i + 1));
  propsy_ = propsy;
  // prepare: zero rules to get constant terms
  lst zeros;
  for (size_t i = 0; i < nint; ++i)
    zeros.append(int_vars_[i] == 0);
  for (size_t i = 0; i < next; ++i)
    zeros.append(ext_vars_[i] == 0);

  // 1. setup the lst sps_to_props_
  //   y = A.x + b
  //   x = A^(-1).(y-b)
  // Steps:
  // (1). compute matrix props_to_sps and constant vector
  // (2). sps_to_props = props_to_sps.inverse()
  // (3). setup the list
  matrix props_to_sps(nprops_, nprops_); // propagators w.r.t. k_i*p_j
  matrix sps_to_props(nprops_, nprops_); // k_i*p_j w.r.t. propagators
  matrix column_vector(nprops_, 1);
  // (1)
  for (size_t i = 0; i < nprops_; ++i) {
    for (size_t j = 0; j < nint; ++j) {
      ex dkj = diff(props_[i], int_vars_[j]);
      // coeff of k_j^2
      props_to_sps(i, j) = diff(dkj, int_vars_[j]) / 2;
      // coeff of k_j*k_l
      size_t prev = nint - 1 + (2 * nint - j - 1) * j / 2;
      for (size_t l = j + 1; l < nint; ++l)
        props_to_sps(i, prev + l - j) = diff(dkj, int_vars_[l]);
      // coeff of k_j*p_l
      for (size_t l = 0; l < next; ++l)
        props_to_sps(i, nint2 + l) = diff(dkj, ext_vars_[l]);
    }
  }
  if (props_to_sps.rank() != nprops_) {
    cout << "\nFailed: not enough irreducible scalar products.\n" << endl;
    exit(-1);
  }
  for (size_t i = 0; i < nprops_; ++i) {
    column_vector(i, 0) = props_[i].subs(zeros, subs_options::algebraic);
  }
  // (2)
  sps_to_props = props_to_sps.inverse();
  // (3)
  column_vector = propsy_.sub(column_vector);
  matrix sps = sps_to_props.mul(column_vector);
  for (size_t j = 0; j < nint; ++j) {
    sps_to_props_.append(int_vars_[j] * int_vars_[j] == sps[j]);
    size_t prev = nint - 1 + (2 * nint - j - 1) * j / 2;
    for (size_t l = j + 1; l < nint; ++l)
      sps_to_props_.append(int_vars_[j] * int_vars_[l] == sps[prev + l - j]);
    for (size_t l = 0; l < next; ++l)
      sps_to_props_.append(int_vars_[j] * ext_vars_[l] == sps[nint2 + l]);
  }
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
  sort(ibps_.begin(), ibps_.end(),
       [](const EquationPtr &l, const EquationPtr &r) {
         return l->integrals.size() < r->integrals.size();
       });

  // save ibp relations to a file
  ofstream ibpfile("./relations/ibp");
  for (size_t i = 0; i < ibps_.size(); ++i)
    ibps_[i].get()->write_file(ibpfile, name_);

  cout << "\nThere are " << ibps_.size() << " IBP relations\n" << endl;
  timer.end();
}

void Tira::create_li_() {
  cout << "\n========== Generate LI relations ===========\n" << endl;
  Timer timer;

  string li_file = "./relations/LI";

  timer.end();
}

void Tira::create_sym_() {}

// IBP relations:
// Integrate[ Diff[k_i^u] * (p_j^u * integrand) ] == 0
// Integrate[ Diff[k_i^u] * (k_j^u * integrand) ] == 0
// Steps:
// 1. general cases: terms from derivatives of the denomimator
// 2. k_j = k_i: an extra term, integrand not changed, coefficient d
// 3. substitute irriducible scalar products and put together
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

      // 1. general cases
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
      // 2. special case when k_j is k_i
      if (vars[j] == int_vars_[i]) {
        RIntegralPtr integral(new RIntegral(nprops_));
        integral->coeff = dimension_;
        ibp->integrals.push_back(move(integral));
      }
      // 3. substitute irriducible scalar products
      subs_irrsp_(ibp);

      // add this relation
      if (ibp->integrals.size() > 0)
        ibps_.push_back(move(ibp));
    }
  }

  delete[] indices;
}

void Tira::subs_irrsp_(EquationPtr &relation) {
  EquationPtr relation_new(new Equation());
  ex coeff;
  sy vv("vv");

  for (auto it = relation->integrals.begin(), ends = relation->integrals.end();
       it != ends; ++it) {
    // substitute irreducible scalar products by propagators
    (*it)->coeff =
        (*it)->coeff.expand().subs(sps_to_props_, subs_options::algebraic);
    (*it)->coeff =
        (*it)->coeff.expand().subs(sp_rules_, subs_options::algebraic);

    // regenerate integrals
    // 1. set dd_k to vv and check if the derivative of vv is zero
    // 2. if not, the integrand is multiplied by prop[k]. Add an new integral
    //    with indices[k]=indices[k]-1 and new coeff;
    // 3. subtract this new coeff from the original coeff (set dd_k to zero)
    // 4. check if the original coeff is zero now
    for (size_t k = 0; k < nprops_; ++k) {
      // 1
      coeff = diff(subs((*it)->coeff.expand(), propsy_(k, 0) == vv,
                        subs_options::algebraic),
                   vv);
      if (coeff != 0) {
        // 2
        RIntegralPtr integral(new RIntegral(nprops_));

        for (size_t j = 0; j < nprops_; ++j)
          integral->indices[j] = (*it)->indices[j];
        --(integral->indices[k]);
        integral->coeff =
            coeff.expand().subs(sp_rules_, subs_options::algebraic);

        relation_new->integrals.push_back(move(integral));
        // 3
        (*it)->coeff =
            (*it)->coeff.subs(propsy_(k, 0) == 0, subs_options::algebraic);
      }
    }
    // 4
    (*it)->coeff = (*it)->coeff.expand();
    if ((*it)->coeff != 0) {
      (*it)->coeff = (*it)->coeff.subs(sp_rules_, subs_options::algebraic);
      relation_new->integrals.push_back(move(*it));
    } else
      (*it).reset();
  }

  // remove duplicates
  relation_new->collect_integrals();

  relation.reset();
  relation = move(relation_new);
}