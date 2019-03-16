#ifndef __Ceiling__
#define __Ceiling__

int Ceiling_log2(int x) {

  unsigned temp = (unsigned) x - 1;
  int result = 0;

  while (temp != 0) {
    temp = temp >> 1;
    result = result + 1;
  }
  return result;
}

#endif
