#include "tira/tools.h"

const possymbol &get_symbol(const string &s) {
  static map<string, possymbol> directory;
  auto i = directory.find(s);
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(make_pair(s, possymbol(s))).first->second;
}