#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include "util.h"

int (* parse(char* input, int* rows))[7] {
  DEBUG_PRINTF("Got input: %s\n", input);
  *rows = -1;
  const int cols = 7;
  char *c = input;

  while (*c != 0) {
    if (*c == '(' || *c == '[' || *c == '{')
      (*rows)++;
    c++;
  }

  DEBUG_PRINTF("Allocating %d rows and %d cols", *rows, 7);

  int (*matrix)[7] = malloc(*rows * 7 * sizeof(int));

  DEBUG_PRINTF("....done\n");

  int row = 0;
  int col = 0;

  for( c = input; *c != 0; c++ ) {

    if (col == 4) {
      /* Currently working under the biological convention*/
      matrix[row][4] = -1;
      /*
       * Un-comment when working under the mathematical convention
       * matrix[row][4] = 1;
      */
      matrix[row][5] = 1;
      matrix[row][6] = 0;
      row++;
      col=0;
      continue;
    }

    if ( *c <= '9' && *c >= '0') {
      matrix[row][col] = strtol(c, &c, 10) - 1;
      col++;
      continue;
    }
  }

#ifdef DEBUG
  printf("DEBUG:\n");
  for (int i = 0; i < *rows; i++){
    for (int j = 0; j < 6; j++){
      printf("%d\t", matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");
#endif

  return matrix;
}
