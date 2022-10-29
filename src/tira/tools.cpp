#include "tira/tools.h"

const possymbol &get_symbol(const string &s) {
  static map<string, possymbol> directory;
  auto i = directory.find(s);
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(make_pair(s, possymbol(s))).first->second;
}

void generate_symbols(sy vars[], const string &s, int n) {
  string str;
  for (unsigned i = 0; i < n; ++i) {
    str = s + to_string(i + 1);
    vars[i] = get_symbol(str);
  }
}