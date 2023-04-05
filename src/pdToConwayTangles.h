#pragma once

void createEdge(int r, int pdCode[r][7], int edge[2 * r][4]);
void getTangle2(int r, int Tangle, int Clock, int *tangle2, int *tangle2Clock,
                int pdCode[r][7], int edge[2 * r][4]);
void combineTanglesEdge(int r, int Tangle, int tangle2, int Arc, int edge[2 * r][4]);
void ClockEdge(int r, int newTang, int Arc, int Clock, int edge[2 * r][4]);
void newArcsCombinedTangle(int r, int pdCode[r][7], int newTangle, int arc0,
                           int arc1, int arc2, int arc3);
int addTangles(int r, int Tangle, int Clock, int Clock2, int pdCode[r][7],
               int edge[2 * r][4]);
int removeTangles(int r, int pdCode[r][7], int edge[2 * r][4], int newRow);
int addRationalTangles(int r, int pdCode[r][7], int edge[][4]);
void rotateTangle(int r, int pdCode[r][7], int tang);
void pdToConway(int r, int pdCode[r][7]);
int compute_writhe(int rows, int pdCode[][7], int edge[2*rows][4] );
