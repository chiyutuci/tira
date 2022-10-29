#include "tira/tira.h"
#include "tira/tools.h"
#include <iostream>

int main(int argc, char **argv) {
  Timer timer;

  Tira tira;
  tira.execute_jobs();

  timer.end();
  return 0;
}