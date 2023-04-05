#pragma once

#ifdef DEBUG
    #define DEBUG_PRINTF(...) printf("DEBUG: "__VA_ARGS__)
#else
    #define DEBUG_PRINTF(...) do {} while (0)
#endif

void display(int row, int cols, int * matrix);

int ainversemodb(int b, int a);

int aModB(int a, int b);
