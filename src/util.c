#include "util.h"
#include <stdio.h>
#include <stdlib.h>

//#ifdef DEBUG
void display(int rows, int cols, int * matrix) {
    for (int i = 0; i < rows; i++){
      for (int j = 0; j < cols; j++){
        printf("%d\t", matrix[i*cols + j]);
      }
      printf("\n");
    }
    printf("\n");
}
/*#else
void display(int row, int cols, int * matrix) {
  do {} while (0);
}
#endif*/
int aModB(int a, int b){
  a = a%b;
  if (a < 0){
    a += b;
  }
  return a;
}
/* ainversemodb finds an X such that X*b = 1 (mod a) */
int ainversemodb(int b, int a) {
  int X = 0;
  if (abs(a)==1){
    return a;
  } else if (a < -1){
    for(int i = -1; i > a; i--){
      if((i * b - 1)  % a == 0){
        X = i - a;
      }
    }
  }
  else {
    for (int i = 1; i < a; i++) {
      if ((i * b - 1) % a == 0){
        X = i;
      }
    }
  }
  if (X == 0){
    printf("ainvermodb broke with values %d, %d\n", b, a);
  }
  return X;
}
