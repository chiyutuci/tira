#ifndef TOOLS_H_
#define TOOLS_H_

#include "ginac/ginac.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <time.h>

using namespace GiNaC;
using namespace std;

typedef possymbol sy;
typedef vector<sy> vsy;

// A simple timer
class Timer {
public:
  Timer() { start = clock(); }

  // time unit: seconds
  double duration() { return (double)(clock() - start) / CLOCKS_PER_SEC; }

  void end() {
    cout << fixed << setprecision(6) << "\nTime: " << duration() << " s\n"
         << endl;
  }
  void end_all() {
    cout << fixed << setprecision(6) << "\nTotal time: " << duration() << " s\n"
         << endl;
  }

private:
  clock_t start;
};

// GiNaC symbol dictionary
const possymbol &get_symbol(const string &s);

// generate a symbol family {s1,s2,...,sn}
void generate_symbols(sy vars[], const string &s, int n);

#endif // TOOLS_H_