#include "tira/tira.h"

void Tira::execute_jobs() {
  config_family();

  create_relations();
}

void Tira::config_family() {
  // register symbols
  // symbols["m"] = get_symbol("m");
  symbols["p1"] = get_symbol("p1");
  symbols["p2"] = get_symbol("p2");
  symbols["p3"] = get_symbol("p3");
  symbols["k1"] = get_symbol("k1");
  parser reader(symbols);

  // family name
  name = "box";

  // dimension d = 4 - 2 * eps
  dimension = get_symbol("d");

  // external momenta
  ext_vars.push_back(get_symbol("p1"));
  ext_vars.push_back(get_symbol("p2"));
  ext_vars.push_back(get_symbol("p3"));

  // internal momenta
  int_vars.push_back(get_symbol("k1"));

  // invariants
  // inv_var=get_symbol("m");

  // invariant set to 1
  // no such symbol

  // scalar products rules
  sp_rules.append(reader("p1") * reader("p1") == reader("0"))
      .append(reader("p1") * reader("p2") == reader("50"))
      .append(reader("p1") * reader("p3") == reader("-1/2"))
      .append(reader("p2") * reader("p2") == reader("0"))
      .append(reader("p2") * reader("p3") == reader("-99/2"))
      .append(reader("p3") * reader("p3") == reader("0"));

  // top level sector, depend on target integrals
  top_sector = 15; // [1,1,1,1]
  // integrand indices, depend on target integrals
  rmax = 6;
  smax = 2;
}

void Tira::create_relations() {}

void Tira::create_ibp() {}

void Tira::create_li() {}

void Tira::create_sym() {}