#include "tira/tira.h"
#include "tira/tools.h"
#include <iostream>

int main(int argc, char **argv) {
  Timer timer;

  Tira tira;
  tira.execute_jobs();

  cout << "\n========== Reduction Finished ==============" << endl;
  timer.end_all();
  return 0;
}