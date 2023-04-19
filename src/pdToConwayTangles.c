#include "util.h"
#include "parse.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pdToConwayTangles.h"
/***************************************************************************************
  Authors: Isabel Darcy, Ethan Rooke, Zachary Bryhtan
  Last Modified: 4-17-2023
  Description: This program takes the PD notation of rational knots then combines the 
    crossings to produce the knot's associated continued fraction. In the 'main' 
    function, the knot's PD notation is stored into a matrix with two additional columns
    used to store a continued fraction for the knot in a given row. A new matrix is 
    created to track which crossings a given arc connects, and where it attaches to 
    each crossing by calling the 'createEdge' function. The 'pdToConway' function is 
    then called which start the process of iterating through the PD notation matrix, 
    'pdCode,' to find crossings which can be combined vertically or horizontally to 
    produce integer tangles. When any two crossings are combined, the associated
    continued fraction is computed to update the values in the last two columns. The 
    row of a merged crossing is filled with zeros in the 'pdCode' matrix. The matrix
    with arc information is also updated to accomodate arcs that are removed after 
    merging, then identifying the new crossings a given arc connects and where it 
    connects. The row of a removed arc is then filled with -1's. After the basic 
    horizontal and vertical combinations are made, the remaining non-zero rows in 
    'pdCode' are iterated through again to find crossings which can be combined, but 
    this time through the resulting crossings will be rational knots. After each 
    combination, the 'pdCode' matrix and the arc matrix are updated to reflect each 
    change. Once there is only one non-zero row in 'pdCode', the last two columns give
    the desired continued fraction. There are several fractions that may be associated
    to a given rational knot, so at this computed fraction is adjusted to provide the
    equivalent fraction with the lowest possible denominator. The continued fraction of
    a knot's mirror image can be readily found from the computed value, so this fraction
    is also given and with the lowest possible denominator.

  Notes: This program also provides an associated sum of coninuted fractions for 
    Montesinos knots and links (hand checked up to 9 crossings). This notation is
    displayed as the numerator closure of the sum of initially found continued fractions
    for each component `N(p1/q1 + ... + pn/qn)', then the sum is transformed into a 
    canonical form as described in a paper by I. Darcy and H. Moon. This form is
    'N(p1/q1 + pn/qn + e)' where now all pi, qi > 0, pi < qi, and e is some integer
    left over after adjusting each summed fraction.
***************************************************************************************/

/* create edge matrix where each ROW corresponds to an edge and
 columns 0 and 2 correspond to the tangles bounding the edge
 columns 1 and 3 indicate if edge corresponds to
 a = 0, b =1, c = 2, d = 3, respectively. */

void createEdge(int r, int pdCode[r][7], int edge[2 * r + 2][4]) {
  for (int i = 0; i < 2 * r + 2; i++){
    for(int j = 0; j < 4; j++){
      edge[i][j] = -9;
    }
  }
  for (int i = 0; i < r; i++){
    for (int j = 0; j < 4; j++) {
      if (edge[pdCode[i][j]][0] == -9) {
        edge[pdCode[i][j]][0] = i;
        edge[pdCode[i][j]][1] = j;
      } else {
        edge[pdCode[i][j]][2] = i;
        edge[pdCode[i][j]][3] = j;
      }
    }
  }
}

void swapTanglesAB(int TangleA, int TangleB, int newRow, int pdCode[][7], int row, int edge[][4]){
  int temp[7];
  for(int i = 0; i < 7; i++){
    temp[i] = pdCode[TangleA][i];
    pdCode[TangleA][i] = pdCode[TangleB][i];
    pdCode[TangleB][i] = temp[i];
  }

  for(int i = 0; i < 2*row; i++){
    if(edge[i][0]==TangleA){
      edge[i][0] = TangleB;
    }else if (edge[i][0]==TangleB){
      edge[i][0] = TangleA;
    }
    if(edge[i][2]==TangleA){
      edge[i][2] = TangleB;
    }else if (edge[i][2]==TangleB){
      edge[i][2] = TangleA;
    }
  }
}

/* When a rational tangle a/b is rotated *by 90 degrees*, it becomes -b/a */
void rotateTangleFraction(int pdCode[][7], int tang) {
  int temp = pdCode[tang][4];
  pdCode[tang][4] = -pdCode[tang][5];
  pdCode[tang][5] = temp;

  if(pdCode[tang][5] < 0){
    pdCode[tang][4] *= -1;
    pdCode[tang][5] *= -1;
  }
}

/*Rotate the given tangle's PD code row and associtated edge information by 90 degrees CCW N-times*/
void rotateTangle90CCW_Ntimes(int N, int Tangle, int totalRows, int pdCode[][7], int edge[][4]){
  int temp[4] = {};
  int rotate = (4 - N)%4;
  for(int i = 0; i < 4; i++){
    temp[i] = pdCode[Tangle][i];
  }
  for(int i = 0; i < 4; i++){
    pdCode[Tangle][i] = temp[(i + rotate) % 4];
  }

  if(rotate % 2 == 1){
    rotateTangleFraction(pdCode, Tangle);
  }

  int edgeRotate = N % 4;
  for(int i = 0; i < 2*totalRows + 2; i++){
    if(edge[i][0] == Tangle){
      edge[i][1] = (edge[i][1] + edgeRotate) % 4;
    }
    if(edge[i][2] == Tangle){
      edge[i][3] = (edge[i][3] + edgeRotate) % 4;
    }
  }

}

/* Given Tangle and Clock, finds *tangle2 which which shares Arc that attaches
 * to Tangle at Clock in {a, b, c, d}.  Also finds *tangle2Clock in {a, b, c, d}
 * where Arc attaches to *tangle2. Using *tangle2 and *tangle2Clock lets us feed
 * initialized variables to the function, update them in the function, then
 * use the result after the call to the function.*/

void getTangle2(int r, int Tangle, int Clock, int *tangle2, int *tangle2Clock,
                int pdCode[][7], int edge[][4]) {
  /* Get the arc label of current tangle through pdCode[Tangle][Clock]
  *  then checks which tangles that arc connect in edge[Arc]
  *  the other tangle index is in the 0, or 2 position followed by the Clock 
  *  position the arc connects through (at position 1 or 3).
  *  We save and change these values in tangle2 and tangle2Clock resp.*/
  int Arc = pdCode[Tangle][Clock];
  if (edge[Arc][0] == Tangle) {
    *tangle2 = edge[Arc][2];
    *tangle2Clock = edge[Arc][3];
  } else if (edge[Arc][2] == Tangle) {
    *tangle2 = edge[Arc][0];
    *tangle2Clock = edge[Arc][1];
  } else
    DEBUG_PRINTF("Tangle = %d, Arc = %d, Arc2 = %d \n", Tangle, edge[Arc][0],
           edge[Arc][2]);
}

//Returns the integer location of another tangle which can be combined to the given tangle
//along the clock positions provided. If no such tangle exists, return -1.
int canCombine(int Tangle, int Clock1, int Clock2, int rows, int pdCode[][7], int edge[][4]){
  int tangle2i, tangle2ii, tangle2iClock, tangle2iiClock;
  
  getTangle2(rows, Tangle, Clock1, &tangle2i, &tangle2iClock, pdCode, edge);
  getTangle2(rows, Tangle, Clock2, &tangle2ii, &tangle2iiClock, pdCode, edge);

  if (tangle2i == tangle2ii && tangle2i >= 0){
    return tangle2i;
  } else {
    return -1;
  }
}

/* changes the edge matrix to denote the new endpoints when tangles are
   combined. If Tangle is eliminated by being combined with tangle2, then
   endpoint Clock of Tangle and its corresponding edge Arc are moved to tangle2.
   Note this subroutine replaces only one endpoint, thus must call twice to
   replace 2  arcs. Clock is replaced using next subroutine, ClockEdge*/

void combineTanglesEdge(int r, int Tangle, int tangle2, int Arc,
                        int edge[2 * r + 2][4]) {
  if (edge[Arc][0] == Tangle) {
    edge[Arc][0] = tangle2;
  }
  if (edge[Arc][2] == Tangle) {
    edge[Arc][2] = tangle2;
  }
}

/* When combining tangles may need to change edge arcs per subroutine below.
   Note Clock is one of (a, b, c, d) = (0, 1, 2, 3).
   Note this might rotate the tangle (handled in addTangles) 
   Updates edge matric to reflect Arc connects to newTang at Clock
*/

void ClockEdge(int r, int newTang, int Arc, int Clock, int edge[2 * r + 2][4]) {
  if (edge[Arc][0] == newTang) {
    edge[Arc][1] = Clock;
  } else if (edge[Arc][2] == newTang) {
    edge[Arc][3] = Clock;
  } else
    DEBUG_PRINTF("mistake");
}

/* places new arcs into pdCode for combined tangle from subroutine addTangles */

void newArcsCombinedTangle(int r, int pdCode[r][7], int newTangle, int arc0,
                           int arc1, int arc2, int arc3) {
  pdCode[newTangle][0] = arc0;
  pdCode[newTangle][1] = arc1;
  pdCode[newTangle][2] = arc2;
  pdCode[newTangle][3] = arc3;
}

//Combines TangleA and TangleB along Clock1 and Clock2 of TangleA, stores result in location of TangleB
//Typically take TangleB < TangleA
int combineTanglesPDandEdge(int TangleA, int TangleB, int Clock1, int Clock2, int rows, int pdCode[][7], int edge[][4]){
  int Clock1B = (Clock1 + 3) % 4;
  int Clock2B = (Clock2 + 1) % 4;
  
  int arc0 = pdCode[TangleA][Clock1B];
  int arc3 = pdCode[TangleA][Clock2B];
  //Swaps TangleA for TangleB in Arc3 row and Arc0 row of Edge
  combineTanglesEdge(rows, TangleA, TangleB, arc0, edge);
  combineTanglesEdge(rows, TangleA, TangleB, arc3, edge);
  
  // ClockEdge(rows, TangleB, arc0, Clock1B, edge);
  // ClockEdge(rows, TangleB, arc3, Clock2B, edge);
  
  pdCode[TangleB][Clock1B] = arc0;
  pdCode[TangleB][Clock2B] = arc3;
  
  //Remove absorbed Arcs from Edge
  for(int i = 0; i < 4; i++){
    edge[pdCode[TangleA][Clock1]][i] = -9;
    edge[pdCode[TangleA][Clock2]][i] = -9;
  }
  //Remove absorbed Tangle from PD
  for(int i = 0; i < 7; i++){
    pdCode[TangleA][i] = 0;
  }
}

/* Determines if a Tangle can be added to the input Tangle.
   Adds these two tangles where the second Tangle is connected
   to the input Tangle via the  two arcs(Arc1, Arc2)
   with endpoints, Clock and Clock2, resp., of the input Tangle.
   returns rotated = 1 iff tangle2i tangle rotated *90* degrees
   Note 180 degree rotation does not change rational tangles,
   but can change non-rational tangles */

int addTangles(int r, int Tangle, int Clock, int Clock2, int pdCode[][7],
               int edge[2 * r+2][4]) {
  
  int temp, j, arc, tangleindex, tangle2i, tangle2ii, tangle2iClock,
      tangle2iiClock=0;
  int rotated = 0;
  int canADD = 0;
  
  // find the index to another tangle connected to this Tangle via the arc at Clock
  tangleindex = Clock;
  getTangle2(r, Tangle, tangleindex, &tangle2i, &tangle2iClock, pdCode, edge);
  
  // find the index to the other tangle connected to this Tangle via the arc at Clock2
  tangleindex = Clock2;
  getTangle2(r, Tangle, tangleindex, &tangle2ii, &tangle2iiClock, pdCode, edge);
  

  // If the above two found tangles are the same, we have found the index to a tangle which shares
  // two arcs with our Tangle, thus should be combinable as a horizontal sum or vertical product
  // There may need to be rotations to complete the task.
  if (tangle2i == tangle2ii && pdCode[tangle2i][4] != 0) {

    /* determine if new tangle2i endpoints require rotation of tangle2i */
    if (tangle2iClock % 2 == 1) {
      rotateTangleFraction(pdCode, tangle2i);
      rotated = 1;
    }

    /* 
    "Add tangles (horizontal sum)"
    PD Code assumes clock positions 0, 1, 2, 3 correspond to SW, SE, NE, NW
    so horizontal sums are really at Clock1 = 3, Clock2 = 0
    or Clock1 = 1 and Clock2 = 2.
     */
      if ((Clock == 3 && Clock2 == 0) || (Clock == 1 && Clock2 == 2)) {
      // addition formula if one of the tangles is (n/1)
      if (abs(pdCode[Tangle][5]) == 1 || abs(pdCode[tangle2i][5]) == 1) {
        canADD = 1;
        pdCode[tangle2i][4] = pdCode[Tangle][4] * pdCode[tangle2i][5] +
                              pdCode[Tangle][5] * pdCode[tangle2i][4];
        pdCode[tangle2i][5] = pdCode[Tangle][5] * pdCode[tangle2i][5];
      }
    }

    /* Multiply tangles (vertical sum) */
      if ((Clock == 3 && Clock2 == 2) || (Clock == 1 && Clock2 == 0)) {
      // addition formula if one of the tangles is (1/n)
      if (abs(pdCode[Tangle][4]) == 1 || abs(pdCode[tangle2i][4] == 1)) {
        canADD = 2;
        pdCode[tangle2i][5] = pdCode[Tangle][4] * pdCode[tangle2i][5] +
                              pdCode[Tangle][5] * pdCode[tangle2i][4];
        pdCode[tangle2i][4] = pdCode[Tangle][4] * pdCode[tangle2i][4];
      }
    }

    if (canADD > 0) {
      /* remove the two interior arcs in the tangle addition */
      if (pdCode[Tangle][Clock] == pdCode[tangle2i][tangle2iClock])
        for (j = 0; j < 4; j++)
          edge[pdCode[Tangle][Clock]][j] = -1;
      else
        DEBUG_PRINTF("mistake1!!, Tangle = %d, tangle2i = %d\n", Tangle, tangle2i);

      if (pdCode[Tangle][Clock2] == pdCode[tangle2ii][tangle2iiClock])
        for (j = 0; j < 4; j++)
          edge[pdCode[Tangle][Clock2]][j] = -1;
      else
        DEBUG_PRINTF("mistake2!! \n");
      /* determine new endpoints/arcs for vertical sum */
      if (Clock == 3 && Clock2 == 2) // Old + New -> New
      {
        /* 
        For matrix pdCode, determine new arcs. The terms 
        (tangle2iClock + 2 )%4 and (tangle2iClock + 3)%4 
        handles rotations of 180 degrees in line.
        */
        newArcsCombinedTangle(r, pdCode, tangle2i,
                              pdCode[Tangle][0],
                              pdCode[Tangle][1],
                              pdCode[tangle2i][(tangle2iClock + 2) % 4],
                              pdCode[tangle2i][(tangle2iClock + 3) % 4]);
        DEBUG_PRINTF("horizontal sum Old + New -> New\n");

        // For matrix edge, determine new arcs coming into new tangle2i.
        // In edge, Tangle is replaced by tangle2i
        combineTanglesEdge(r, Tangle, tangle2i, pdCode[Tangle][0], edge);
        combineTanglesEdge(r, Tangle, tangle2i, pdCode[Tangle][1], edge);
        // In edge, clock of tangle2i is changed
        ClockEdge(r, tangle2i, pdCode[tangle2i][2], 2, edge);
        ClockEdge(r, tangle2i, pdCode[tangle2i][3], 3, edge);
      }

      /* determine new endpoints/arcs for vertical sum*/
      if (Clock == 1 && Clock2 == 0) // New + Old -> New
      {
        // For matrix pdCode, determine new arcs.
        newArcsCombinedTangle(r, pdCode, tangle2i,
                              pdCode[tangle2i][(tangle2iClock + 2) % 4],
                              pdCode[tangle2i][(tangle2iClock + 3) % 4],
                              pdCode[Tangle][2], pdCode[Tangle][3]);
        DEBUG_PRINTF("horizontal sum2 New + Old -> New \n");

        // For matrix edge, determine new arcs coming into new tangle2i.
        // In edge, Tangle is replaced by tangle2i
        combineTanglesEdge(r, Tangle, tangle2i, pdCode[Tangle][2], edge);
        combineTanglesEdge(r, Tangle, tangle2i, pdCode[Tangle][3], edge);
        // In edge, clock of tangle2i is changed
        ClockEdge(r, tangle2i, pdCode[tangle2i][0], 0, edge);
        ClockEdge(r, tangle2i, pdCode[tangle2i][1], 1, edge);
      }

      /* determine new endpoints/arcs for horizontal sum */
      if (Clock == 1 && Clock2 == 2) // Old * New -> New
      {
        newArcsCombinedTangle(r, pdCode, tangle2i, pdCode[Tangle][0],
                              pdCode[tangle2i][(tangle2iClock + 1) % 4],
                              pdCode[tangle2i][(tangle2iClock + 2) % 4],
                              pdCode[Tangle][3]);
        DEBUG_PRINTF("vertical sum Old * New -> New \n");

        // Determine new arcs coming into new tangle2i in matrix edge
        // In edge, Tangle is replaced by tangle2i
        combineTanglesEdge(r, Tangle, tangle2i, pdCode[Tangle][0], edge);
        combineTanglesEdge(r, Tangle, tangle2i, pdCode[Tangle][3], edge);
        // In edge, clock of tangle2i is changed
        ClockEdge(r, tangle2i, pdCode[tangle2i][1], 1, edge);
        ClockEdge(r, tangle2i, pdCode[tangle2i][2], 2, edge);
      }

      /* determine new endpoints/arcs for horizontal sum */
      if (Clock == 3 && Clock2 == 0) // New*Old -> New
      {
        newArcsCombinedTangle(r, pdCode, tangle2i,
                              pdCode[tangle2i][(tangle2iClock + 2) % 4],
                              pdCode[Tangle][1], pdCode[Tangle][2],
                              pdCode[tangle2i][(tangle2iClock + 1) % 4]);
        DEBUG_PRINTF("vertical sum ");

        // Determine new arcs coming into new tangle2i in matrix edge
        // In edge, Tangle is replaced by tangle2i
        combineTanglesEdge(r, Tangle, tangle2i, pdCode[Tangle][1], edge);
        combineTanglesEdge(r, Tangle, tangle2i, pdCode[Tangle][2], edge);
        // In edge, clock of tangle2i is changed
        ClockEdge(r, tangle2i, pdCode[tangle2i][0], 0, edge);
        ClockEdge(r, tangle2i, pdCode[tangle2i][3], 3, edge);
      }

      // remove Tangle
      for (j = 0; j < 7; j++)
        pdCode[Tangle][j] = 0;

      DEBUG_PRINTF("Tangle = %d \n", Tangle);

      DEBUG_PRINTF("Clock = %d, tangle2i = %d, and its clock = %d \n", Clock,
             tangle2i, tangle2iClock);
    }
  }
  if (rotated == 1 && canADD == 0) {
    rotateTangleFraction(pdCode, tangle2i);
    rotated = 1;
    pdCode[tangle2i][4] = -pdCode[tangle2i][4];
    pdCode[tangle2i][5] = -pdCode[tangle2i][5];
    DEBUG_PRINTF("NOT rotated\n");
  }
  return canADD;
  //return rotated;
}

// remove rows in pdCode where tangle has been combined with another tangle
/*
  Slides rows in pdCode up to move all zero rows to the bottom. The edge matrix
  is also updated to accomodate these changes
*/
int removeTangles(int r, int pdCode[r][7], int edge[2 * r + 2][4], int newRow) {
  int i, j, k;
  //lowers each crossing label in edge for each arc by the number of zero rows in pdCode
  // for loop must be descending to ensure we reduce the crossing number correctly
  for (i = newRow; i > 0; i--) {
    if (pdCode[i][4] == 0 || pdCode[i][5] == 0) {
      for (j = 0; j < 2 * r + 2; j++) {
        if (edge[j][0] > i - 1){
          edge[j][0] = edge[j][0] - 1;
        }
        if (edge[j][2] > i - 1){
          edge[j][2] = edge[j][2] - 1;
        }
      }
    }
  }

  int Start = newRow;
  for (i = Start - 1; i > -1; i--) {
    if (pdCode[i][4] == 0 || pdCode[i][5] == 0) {
      if (i == newRow - 1)
        newRow = newRow - 1;
      else {
        for (j = i; j < newRow - 1; j++)
          for (k = 0; k < 7; k++) {
            pdCode[j][k] = pdCode[j + 1][k];
          }
        for (k = 0; k < 7; k++) {
          pdCode[newRow - 1][k] = 0;
        }
        newRow = newRow - 1;
      }
    }
  }
  return newRow;
}


/* 
If newRow = 2 and pdCode contains 2 rational tangles
Then tangle is rational since
N(a/b + t/w) = N((aw + bt)/(xw + yt))
             = N((aw + tb)/(tb' + a'w))
             = N((aw + tb)/(ty  + xw))
where a'b - ab' = 1 and x b - ay = 1
Thus a' = x and b' = y
*/
int addRationalTangles(int r, int pdCode[r][7], int edge[][4]) {
  int a, b, t, w, x, y;
  
  //Store current numerator and denominator of first tangle as 'a' and 'b' respectively.
  a = pdCode[0][4];
  b = pdCode[0][5];
  // Store current numerater and denominator of second tangle as 't' and 'w' respectively.
  t = pdCode[1][4];
  w = pdCode[1][5];


  if (pdCode[0][0] == pdCode[1][3] && pdCode[0][1] == pdCode[1][2] &&
      pdCode[0][2] == pdCode[1][1] && pdCode[0][3] == pdCode[1][0]) {
    // Horizontal sum could work either T1 + T2 or T2 + T1
    /* 
    Removing these lines as clutter, but this whole conditional might be unnecessary
    a = pdCode[0][4];
    b = pdCode[0][5];
    t = pdCode[1][4];
    w = pdCode[1][5];
    */
  } else if (pdCode[0][0] == pdCode[1][1] && pdCode[0][1] == pdCode[1][0] &&
             pdCode[0][2] == pdCode[1][3] && pdCode[0][3] == pdCode[1][2]) {
    // Vertical product could work either T1 * T2 or T2 * T1
    /*
    Indicating a rotation by 90 degrees of both tangles in their fraction by taking
    each n/m to m/(-n).
    */
    int temp = b;
    b = -a; // vertical sum
    a = temp;  // Note a/b rotated 90 degrees = -b/a
    temp = w;
    w = -t;
    t = temp;
  } else {
    //make sure tangle containing 0 is first
    if(edge[0][0] != 0){
      swapTanglesAB(0, 1, 2, pdCode, r, edge);
    }

    if(edge[0][1]!= 0){
      rotateTangle90CCW_Ntimes(4 - edge[0][1], 0, r, pdCode, edge);
    }
      //Now the tangle containing '0' edge is oriented with 0 in SW
    
    if(edge[pdCode[0][1]][2] < 0){
      rotateTangle90CCW_Ntimes(3, 0, r, pdCode, edge);
    }

    // if(pdCode[0][5] < 0){
    //   pdCode[0][4] *= -1;
    //   pdCode[0][5] *= -1;
    // }

    a = pdCode[0][4];
    b = pdCode[0][5];
   
    //bot_left tangle is positioned to be summed with the remaining component
    //rotate top_right so that the edge label at pdCode[bot_left][1] is at
    //clock 0 in the top_right tangle

    int arc1 = pdCode[0][1];
    if( edge[arc1][0] == 1 && edge[arc1][1] != 0){
      rotateTangle90CCW_Ntimes(4 - edge[arc1][1], 1, r, pdCode, edge);
    } else if (edge[arc1][2] == 1 && edge[arc1][3] != 0){
      rotateTangle90CCW_Ntimes(4 - edge[arc1][3], 1, r, pdCode, edge);
    }

    // if(pdCode[top_right][5] < 0){
    //   pdCode[top_right][4] *= -1;
    //   pdCode[top_right][5] *= -1;
    // }

    t = pdCode[1][4];
    w = pdCode[1][5];    

  }
  
  // the setup is a sum, if the second term is integral, combination is easy
  if( w == 1){
    pdCode[0][4] = a + b*t;
    pdCode[0][1] = pdCode[1][1];
    pdCode[0][2] = pdCode[1][2];
    
    for(int i = 0; i < 4; i++){
      edge[pdCode[1][0]][i] = -9;
      edge[pdCode[1][3]][i] = -9;
    }
    for(int i = 0; i < 7; i++){
      pdCode[1][i] = 0;
    }

    return 1;
  }
  

  x = ainversemodb(b, a);
  y = (b * x - 1) / a;
  /* 
  This check is probably uneeded
  if (b * x - a * y != 1)
    DEBUG_PRINTF("mistake");
  */
  /*
  ''correct this, should be closure with arcs''
  Not sure what this means, but the loop is clearly changing nothing
  for (int i = 0; i < 4; i++){
    pdCode[0][i] = pdCode[0][i];
  }
  */
  pdCode[0][4] = a * w + b * t;
  pdCode[0][5] = x * w + y * t;
  if (abs(pdCode[0][4]) < abs(pdCode[0][5])){
    pdCode[0][5] = pdCode[0][5]%(abs(pdCode[0][4]));
  }
  int oneRow = 1;
  return oneRow;
}

void insertZeroRows(int where, int amount, int newRow, int pdCode[][7]){
  if(amount != 0){
    for(int i = newRow + amount; i > where + amount; i--){
      for(int j = 0; j < 7; j++){
        pdCode[i][j] = pdCode[i - amount][j];
      }
    }
    for(int i = where+1; i < where + amount; i++){
      for(int j = 0; j < 7; j++){
        pdCode[i][j] = 0;
      }
    }
  }
}

/*char combineComponents(int rows, int currentComponent, int *component2 int reducedPdCode[][4], edge[][4]){
  int tangle2i,tangle2ii,tangle2iClock, tangle2iiClock = 0;
  getTangle2(rows, currentComponent, 0, &tangle2i, &tangle2iClock, reducedPdCode, edge);
  getTangle2(rows, currentComponent, 1, &tangle2ii, &tangle2iiClock, reducedPdCode, edge);
  if(tangle2i == tangle2ii){
    if(tangle2iClock % 2 == 1){
      int adjust = (tangle2iiClock + 2)%4;
      for(int i = 0; )
    }
    *component2 = tangle2i;
    return '*';
  }
}*/

int main(int argc, char *argv[]) {
  if ( argc < 2) return 1;
  
  int row;
  int (*pdCode)[7] = parse(argv[1], &row);
  //The seventh column stores operation so that
  //tangle i connects to tangle i+1 by operation pdCode[i][6]
  //operations are -1 = *, 0 = ?, 1 = +.

  pdToConway(row, pdCode);
}

/*
 * Given the pdCode of a knot along with its corresponding edge matrix
 * compute the writhe of the knot.
 */
int compute_writhe(int rows, int pdCode[][7], int edge[2*rows + 2][4] ){

  int crossing2 = 0;
  int crossing2Clock = 0;
  int writhe = 0;
  /*Assuming the PD notation given follows expected labeling conventions, we need 
   *only take the difference of edge labels (Clock 3) - (Clock 1). Otherwise, we
   *may need the following procedure, but its debatable.*/

  for (int i = 0; i < rows; i++) {
    /* We definitely don't need both calls searching clocks 1 AND 3, just one is
     * sufficient. In fact the search at clock 3 in general overwrites the initial
     * search with exceptions for an edge which is an understrand through two 
     * crossings*/
    int sign = 0;
    /* Looking for second tangle (crossing2) which connects to tangle 'i' via the 
    * edge at its clock position 1 (bottom right) and at which position (crossing2Clock).*/
    getTangle2(rows, i, 1, &crossing2, &crossing2Clock, pdCode, edge);

    if (crossing2Clock == 0){
      sign = -1;
    }
    else if (crossing2Clock == 2){
      sign = 1;
    }

    /* Looking for second tangle (crossing2) which connects to tangle 'i' via the
     * edge at clock position 3 (top left) and at which position (crossing2Clock) */
    getTangle2(rows, i, 3, &crossing2, &crossing2Clock, pdCode, edge);

    if (crossing2Clock == 0){
      sign = 1;
    }
    else if (crossing2Clock == 2){
      sign = -1;
    }

    if (sign == 0)
      DEBUG_PRINTF("incomplete");

    writhe += sign;
  }

  return writhe;
}

void sort(int row, int newRow, int pdCode[][7], int edge[][4]){
    
  int temp[6];
  int sorted = 1;
  while(sorted > 0){
    for (int i=0; i < newRow-1; i++){
      sorted = 0;
      if(abs(pdCode[i][4]) < abs(pdCode[i+1][4])){
        sorted = 1;
        for(int j = 0; j < 6; j++){
          temp[j] = pdCode[i][j];
          pdCode[i][j] = pdCode[i+1][j];
          pdCode[i+1][j] = temp[j];
        }
        for(int k = 0; k < 2*row; k++){
          if(edge[k][0] == i){
            edge[k][0] = i+1;
          } else if (edge[k][0] == i+1){
            edge[k][0] = i;
          }
          if(edge[k][2] == i){
            edge[k][2] = i+1;
          }else if (edge[k][2] == i+1){
            edge[k][2] = i;
          }
        } 
      }
    }
  }
}

void makeCanonical(int start, int end, int pdCode[][7], int sign, int *remainder, int orientation){

  if(end - start > 1){
    if(sign == -1){
      printf("-");
      for (int i = start; i < end; i++){
        pdCode[i][4] *= -1;
      }
      *remainder *= -1;
    }
    if(orientation%2==0){
      int temp = 0;
      for(int i = start; i < end; i++){
        temp = pdCode[i][4];
        pdCode[i][4] = pdCode[i][5];
        pdCode[i][5] = temp;
      }
      //*remainder *= -1;
      
    }
    
    
    int a, b, q;  
    for(int i = start; i < end; i++){
      a = pdCode[i][4];
      b = pdCode[i][5];
      if(b == 1 || a == 0){
        /*mostly unecessary check assuming all integer/vertical/rational tangles
        *were created correctly, and the isMontesinos check passed.*/
        printf("Potentially non-Montesinos");
        exit(-1);
      } else {
        q = aModB(a, b);
        pdCode[i][4] = q;
        *remainder += (a - q)/b;
      }
    }
    
    // We should sort by largest separation (den) - (num)
    int temp = 0;
    for(int i = start; i < end-1; i++){
      if(pdCode[i][5] - pdCode[i][4] < pdCode[i+1][5] - pdCode[i+1][4]||pdCode[i][5] < pdCode[i+1][5]){
        for(int j = 0; j < 6; j++){
          temp = pdCode[i][j];
          pdCode[i][j] = pdCode[i+1][j];
          pdCode[i+1][j] = temp;
        }
      }
    }
    if(orientation%2==0){
      int temp = 0;
      for(int i = start; i < end; i++){
        temp = pdCode[i][4];
        pdCode[i][4] = pdCode[i][5];
        pdCode[i][5] = temp;
      }
    }
  } else if (end - start == 1&&sign == -1){
    printf("-");
    pdCode[start][4] *= -1;

  }
}
/*
  Accepts pdCode and operations to check whether the algebraic tangle
  a Montesinos tangle, i.e. may be expressed as either a horizontal
  sum of vertical tangles or as a verticle product of horizontal tangles
*/
int isMontesino(int startRow, int endRow, int pdCode[][7], int isComponent){
  //-1 is always bad, we need an operation unless its polyhedral
  if(pdCode[startRow][6] == 0){
    return -1;
  }
  //plan to use the +1 value to flag an algebraic tangle
  if(pdCode[startRow][6] == 1){
    for(int i = startRow + 1; i < endRow - 1; i++){
      if(pdCode[i][6] == -1||pdCode[i][6] == 0){
        return -1;
      }
    }
    return 0;
  }
  if(pdCode[startRow][6] == -1){
    for(int i = startRow + 1; i < endRow - 1; i++){
      if(pdCode[i][6] == 1||pdCode[i][6] == 0){
        return -1;
      }
    }
    if(isComponent == 0){
      /*escaping the loop means we identified only products, so we go
      *through and rotate the diagram by 90 degrees to match expected
      *sums for Montesinos structures.*/

      for (int i = startRow; i < endRow; i++){
        rotateTangleFraction(pdCode, i);
        if(i < endRow-1){
          pdCode[i][6] = 1;
        }
      }
      //clean up negative signs
      for(int i = startRow; i < endRow; i++){
        if(pdCode[i][5] < 0){
          pdCode[i][4] *= -1;
          pdCode[i][5] *= -1;
        }
      }
      return 0;
    }
    return 1;
  }
}

int majoritySign(int start, int end, int pdCode[][7], int *remainder, int orientation){
  int posCount=0;
  int adjust=0;
  if(orientation%2==0){
    for(int i = start; i < end; i++){
      //pdCode[i][4]*= -1;
      rotateTangleFraction(pdCode, i);
      if(pdCode[i][4] < 0 && pdCode[i][5] < 0){
        pdCode[i][4] *= -1;
        pdCode[i][5] *= -1;
      }
    }
  }
  for(int j = start; j < end; j++){
    if(pdCode[j][4]>0){
      posCount++;
    }
  }
  if(posCount > (end - start)/2){
    //if majority is +, identify - and adjust
    for(int j = start; j < end; j++){
      if(pdCode[j][4] < 0){
        adjust = -pdCode[j][4]/pdCode[j][5] + 1;
        *remainder -= adjust;
        pdCode[j][4] += adjust*pdCode[j][5];
      }
    }
    if(orientation%2==0){
      for(int i = start; i < end; i++){
        rotateTangleFraction(pdCode, i);
        if(pdCode[i][4] < 0 && pdCode[i][5] < 0){
          pdCode[i][4] *= -1;
          pdCode[i][5] *= -1;
        }
      }
      *remainder *= -1;
      return -1;
    }
    return 1;
  } else {
    
    //if majority is -, identify + and adjust
    for(int j = start; j < end; j++){
      if(pdCode[j][4] > 0){
        adjust = -pdCode[j][4]/pdCode[j][5] - 1;
        *remainder -= adjust;
        pdCode[j][4] +=adjust*pdCode[j][5];
      }
    }
    if(orientation%2==0){
      for(int i = start; i < end; i++){
        rotateTangleFraction(pdCode, i);
        if(pdCode[i][4] < 0 && pdCode[i][5] < 0){
          pdCode[i][4] *= -1;
          pdCode[i][5] *= -1;
        }
      }
      
      *remainder *=-1;
      return 1;
    }
    return -1;
  }
}

void getFraction(int start, int end, int newRow, int pdCode[][7], int *remainder, int sign, int orientation){
  if(end - start == newRow){
    printf("N(%d/%d", pdCode[0][4], pdCode[0][5]);
    for(int i = 1; i < newRow; i++){
      if(pdCode[i-1][6]==1){
        printf(" + ");
      }
      if(pdCode[i-1][6]==-1){
        printf(" * ");
      }
      printf("%d/%d", pdCode[i][4], pdCode[i][5]);
    }
    if(*remainder!= 0){
      printf(" + %d),", *remainder);
    } else {
      printf("),");
    }
  } else{
    printf("(%d/%d", pdCode[start][4], pdCode[start][5]);
  
    for(int i = start; i < end-1; i++){
      if(pdCode[i][6]==1){
        printf(" + ");
      }
      if(pdCode[i][6]==-1){
        printf(" * ");
      }
      printf("%d/%d", pdCode[i+1][4], pdCode[i+1][5]);
    
    }
  
    if(*remainder == 0){
      printf(")");
    } else {
      if(orientation%2==1){
        printf(" + %d)", *remainder);
      }
      else{
        printf(" * 1/%d)", *remainder);
      }
    }
  }
}

void getConway(int a, int b){
  int A_i = a/b; //Conway output this step
  int r = aModB(a, b); // a/b = A_i + r/b
  if(r > 1){//if remainder is more than one, flip the fraction and do it again
    getConway(b, r);
  } else if (r > 0){
    printf("%d", b);
  }
  printf(" %d", A_i);
}

void getConwayMontesinos(int sign,int remainder, int start, int stop, int newRow, int pdCode[][7]){
  
  if(sign == -1){
        printf("-");
      }
      if(stop - start < newRow){
        printf("(");
      }else {
        printf("[");
      }
      for(int i = start; i < stop; i++){
        getConway(pdCode[i][4], pdCode[i][5]);
        if(i < stop - 1){
          printf("|");
        }
      }
      if(remainder != 0){
        char pm = '+';
        if(remainder < 0){
          pm = '-';
        }
        for(int i = 0; i < abs(remainder); i++){
          printf("%c", pm);
        }
      }
      if(stop - start < newRow){
        printf(")");
      } else{
        printf("]");
      }
}

void handleMontesinosComponent(int start, int end, int pdCode[][7], int mont){
  printf("(%d/%d", pdCode[start][4], pdCode[start][5]);
  for(int i = start+1; i < end; i++){
    if(mont == 0){
      printf(" + ");
    } else if (mont == 1){
      printf(" * ");
    }
    printf("%d/%d", pdCode[i][4], pdCode[i][5]);
  }
  printf(")");
}

void handleMontesinos(int newRow, int pdCode[][7]){
  int mont = isMontesino(0, newRow, pdCode, 0);
    if(mont == 0){
      //Assuming Montesinos, we output current continued fraction expression
      printf("N(%d/%d", pdCode[0][4], pdCode[0][5]);
      for(int i = 0; i < newRow-1; i++){
        if(pdCode[i+1][4]!=0){
          if(pdCode[i][6] == 1){
            printf(" + ");
          }
          if(pdCode[i][6]==-1){
            printf(" * ");
          }
          printf("%d/%d", pdCode[i+1][4], pdCode[i+1][5]);
        }
      }
      printf("),");
      //Now we want to get each fraction into the same sign, whatever is the majority
      int remainder=0;
      int sign = 0;
      sign = majoritySign(0, newRow, pdCode, &remainder, 1);
      makeCanonical(0, newRow, pdCode, sign, &remainder, 1);
      getFraction(0, newRow,newRow, pdCode, &remainder, sign, 1);
      getConwayMontesinos(sign, remainder,0, newRow, newRow, pdCode);
      printf(",");

      getFraction(0, newRow, newRow, pdCode, &remainder, -sign, 1);
      getConwayMontesinos(-sign, remainder, 0, newRow, newRow, pdCode);
      printf(",");
    }
}

void orientAlgebraic(int pieces, int newRow, int row, int components[], int pdCode[][7], int edge[][4]){
  // if pieces is even, first piece should be oppositely oriented
      //i.e. expected mont operation is * or (-1)
    // if pieces is odd, first piece should be normally oriented
      //i.e. expected mont operation is + or (1)
  int tangle2i, tangle2ii, tangle2iClock, tangle2iiClock;
  int temp[7] = {0,0,0,0,0,0,0};
  int top_right = 0; //component where top-right tangle lives
  int bottom_left = 0; //component where bottom-left tangle lives
    /* int r1[7];
  const int bytes = 7*sizeof(int);
  memcpy(r1, pdCode[0], bytes);
  memcpy(pdCode[0], pdCode[1], bytes);
  memcpy(pdCode[1], r1, bytes); */
  
  //Now start by looking for a piece that combines with the first
  int previous_connection = 0;
  for(int i = 0; i < pieces-1; i++){

    if(i > 0){
      previous_connection = pdCode[components[i]-1][6];
    }

    if((i==0 && pdCode[0][6]!=-1) || previous_connection==1){//then this is a summed component
      //as a summed component, we should look for a product first
      //product with current stack useslast tangle
      
      int left_tangle = components[bottom_left];
      int right_tangle = components[top_right+1]-1;
      if(pdCode[components[bottom_left]][6]==-1){
        left_tangle = components[bottom_left+1]-1;
      }
      getTangle2(row, left_tangle, 0, &tangle2i, &tangle2iClock, pdCode, edge);
      getTangle2(row, right_tangle, 1, &tangle2ii, &tangle2iiClock, pdCode, edge);
      
      int here = 0;
     
      for(int j = i; j < pieces; j++){
        if(tangle2i < components[j+1] && tangle2ii < components[j+1]){
          here = j;
          break;
        }
      }
      if(here!=0&&((tangle2i==tangle2ii)||(tangle2i==components[here] && tangle2ii==components[here+1]-1)||(tangle2ii==components[here] && tangle2i==components[here + 1] - 1))){
        //Then the tangle(s) belong to the same component so we can proceed
        
        int rotate = tangle2iiClock;//2 is none, 3 is 90 cw, 0 is 180, 1 is 90 ccw
        if(here > i+1){//Then we need to swap these components
          int size = 0;
          int diff = 0;
          if(components[here+1] - components[here] >= components[i+2] - components[i+1]){
            size = components[i+2] - components[i+1];
            diff = components[here+1]-components[here] - size;
            insertZeroRows(components[i+2]-1, diff, newRow, pdCode);
            //included zero rows to components[i+1] piece to make them equal sizes
            for(int j = components[i+1]; j < components[i+1] + size + diff; j++){
              swapTanglesAB(j, components[here] + j - components[i+1], newRow, pdCode,row, edge);
            }
            for(int j = components[here+1]-diff;j < newRow+diff; j++){
              for(int k = 0; k < 7; k++){
                pdCode[j][k] = pdCode[j+diff][k];
              }
            }

          } else {
            size = components[here+1] - components[here];
            diff = components[i+2] - components[i+1] - size;
            insertZeroRows(components[here+1]-1, diff, newRow, pdCode);
            //included zero rows to components[location] piece to make them equal sizes
            for(int j = components[i+1]; j < components[i+1] + size + diff; j++){
              swapTanglesAB(j, components[here] + j - components[i+1], newRow, pdCode,row, edge);
            }
            
            for(int j = components[i+2]-diff; j < newRow + diff; j++){
              for(int k = 0; k < 7; k++){
                pdCode[j][k] = pdCode[j+diff][k];
                if((k < 4) && edge[pdCode[j][k]][0]==j+diff){
                  edge[pdCode[j][k]][0] -= diff;
                }
                if((k < 4) && edge[pdCode[j][k]][2]==j+diff){
                  edge[pdCode[j][k]][2] -= diff;
                }
              }
            }
          }
          
          

        //Components are swapped, one at location is now at i+1*******************
          //now update components array
          
          int x = 1;
          for(int j = 1; j < newRow; j++){
            if(pdCode[j-1][6]==0){
              components[x] = j;
              x++;
            }
          }
          
        }
        
        
        //Components are swapped and component array is updated
        //rotate component as needed to perform the operation
        if(rotate%2==1){//odd rotations need flipped fractions a/b -> -b/a
          for(int j = components[i+1]; j < components[i+2]; j++){
            rotateTangleFraction(pdCode, j);
            for(int k = 0; k < 4; k++){
              temp[k] = pdCode[j][k];
            }
            for(int k = 0; k < 4; k++){
              //rotate = 1 is 90 CCW (0->1->2->3->0); orig -> (orig + rotate)%4
              //rotate = 3 is 90 CW (0->3->2->1->0); orig -> (orig + rotate)%4
              pdCode[j][k] = temp[(k + 2 + rotate)%4];
            }
            for(int k = 0; k < 2*row; k++){
              if(edge[k][0]==j){
                edge[k][1] = (edge[k][1] + 2 + rotate)%4;
              }
              if(edge[k][2]==j){
                edge[k][3] = (edge[k][3] + 2 + rotate)%4;
              }
            }
            pdCode[j][6]*=-1;
          }
          
          if((rotate==1 && pdCode[components[i+1]][6] == -1)||(rotate == 3 && pdCode[components[i+1]][6] == 1)){//CCW makes first in sum, last in product
            temp[0] = pdCode[components[i+1]][6];
            pdCode[components[i+1]][6] = pdCode[components[i+2]-1][6];
            pdCode[components[i+2]-1][6] = temp[0];
            int k = 1;
            for(int j = components[i+1]; j < (components[i+2]+components[i+1])/2; j++){
              swapTanglesAB(j, components[i+2]-k, newRow, pdCode, row, edge);
              k++;
            }
          }
        }
        else if (rotate == 0){
          temp[0] = pdCode[components[i+1]][6];
          pdCode[components[i+1]][6] = pdCode[components[i+2]-1][6];
          pdCode[components[i+2]-1][6] = temp[0];
          int k = 1;
          for(int j = components[i+1]; j < (components[i+2]+components[i+1])/2; j++){
            swapTanglesAB(j, components[i+2]-k, newRow, pdCode, row, edge);
            k++;
          }
          for(int j = components[i+1]; j < components[i+2]; j++){
            for(int k = 0; k < 4; k++){
              temp[k] = pdCode[j][k];
            }
            for(int k = 0; k < 4; k++){
              pdCode[j][k] = temp[(k + 2)%4];
            }
            for(int k = 0; k < 2*row; k++){
              if(edge[k][0]==j){
                edge[k][1] = (edge[k][1] + 2)%4;
              }
              if(edge[k][2]==j){
                edge[k][3] = (edge[k][3] + 2)%4;
              }
            }
          }
        }

        pdCode[components[i+1]-1][6] = -1;
        pdCode[components[i+2] - 1][6] = 0;
        bottom_left = i+1;
        continue;
      }

      //sum with component is at last tangle, 2nd clock and last tangle 1st clock
      int top_tangle = components[top_right];
      int bottom_tangle = components[top_right + 1] - 1;
      if(pdCode[components[top_right]][6]==1){
        top_tangle = bottom_tangle;
      }
      getTangle2(row, top_tangle, 2, &tangle2i, &tangle2iClock, pdCode, edge);
      getTangle2(row, bottom_tangle, 1, &tangle2ii, &tangle2iiClock, pdCode, edge);
      
      here = 0;
      for(int j = i; j < pieces; j++){
        if(tangle2i < components[j+1] && tangle2ii < components[j+1]){
          here = j;
          break;
        }
      }
      //We have found which component contains tangle2i and tangle2ii
      //We now have to decide whether the component needs to be moved so that it occurs next in pdCode
      //Then decide whether the component must be rotated
      if(here != 0&&((tangle2i==tangle2ii)||(tangle2i==components[here] && tangle2ii==components[here + 1]-1)||(tangle2ii==components[here]&&tangle2i==components[here + 1]-1))){
        //Then the tangles belong to the same component to we can proceed
        int rotate = tangle2iiClock;//0 is none, 1 is 90 cw, 2 is 180, 3 is 90 ccw
        if(here > i+1){//Then we need to swap these components        
          int size = 0;
          int diff = 0;
          if(components[here+1] - components[here] >= components[i+2] - components[i+1]){
            size = components[i+2] - components[i+1];
            diff = components[here+1]-components[here] - size;
            insertZeroRows(components[i+2]-1, diff, newRow, pdCode);
            //included zero rows to components[i+1] piece to make them equal sizes
            for(int j = components[i+1]; j < components[i+1] + size + diff; j++){
              swapTanglesAB(j, components[here] + j - components[i+1], newRow, pdCode,row, edge);
            }
            for(int j = components[here+1]-diff;j < newRow+diff; j++){
              for(int k = 0; k < 7; k++){
                pdCode[j][k] = pdCode[j+diff][k];
              }
            }
          } else {
            size = components[here+1] - components[here];
            diff = components[i+2] - components[i+1] - size;
            insertZeroRows(components[here+1]-1, diff, newRow, pdCode);
            //included zero rows to components[location] piece to make them equal sizes
            for(int j = components[i+1]; i < components[i+1] + size + diff; j++){
              swapTanglesAB(j, components[here] + j - components[i+1], newRow, pdCode,row, edge);
            } 
            
            for(int j = components[i+2]-diff; j < newRow + diff; j++){
              for(int k = 0; k < 7; k++){
                pdCode[j][k] = pdCode[j+diff][k];
                if((k < 4) && edge[pdCode[j][k]][0]==j+diff){
                  edge[pdCode[j][k]][0] -= diff;
                }
                if((k < 4) && edge[pdCode[j][k]][2]==j+diff){
                  edge[pdCode[j][k]][2] -= diff;
                }
              }
            }
          }
            //Components are swapped, one at location is now at i+1*******************
            //now update components array
          int x = 1;
          for(int j = 1; j < newRow; j++){
            if(pdCode[j-1][6]==0){
              components[x] = j;
              x++;
            }
          }
        }
      
        //****Now rotate component in pdCode as needed***************************
        if(rotate%2 == 1){//odd rotations need flipped fractions a/b -> -b/a
          for(int j = components[i+1]; j < components[i+2] - 1; j++){
            rotateTangleFraction(pdCode, j);
            for(int k = 0; k < 4; k++){
              temp[k] = pdCode[j][k];
            }
            for(int k = 0; k < 4; k++){
              //rotate = 1 is 90 CW (0->3->2->1->0); orig -> (orig + 2 + rotate) mod 4
              //rotate = 3 is 90 CCW (0->1->2->3->0); orig -> (orig + 2 + rotate) mod 4
              pdCode[j][k] = temp[(k + rotate)%4];
            }
            for(int k = 0; k < 2*row; k++){
              if(edge[k][0]==j){
                edge[k][1] = (edge[k][1] +rotate)%4;
              }
              if(edge[k][2]==j){
                edge[k][3] = (edge[k][3] +rotate)%4;
              }
            }
            pdCode[j][6] *= -1;
          }
          if ((rotate==1&&pdCode[components[i+1]][6]==1)||(rotate == 3 && pdCode[components[i+1]][6] == -1)){//reverse the order of tangles making up the Montesinos component
            temp[0] = pdCode[components[i+1]][6];
            pdCode[components[i+1]][6] = pdCode[components[i+2]-1][6];
            pdCode[components[i+2]-1][6] = temp[0];
            int k = 1;
            for(int j = components[i+1]; j < (components[i+2]+components[i+1])/2; j++){
              swapTanglesAB(j, components[i+2]-k, newRow, pdCode,row, edge);
              k++;
            }
          }
        }
        else if (rotate==2){//could also maybe do rotate%2==0, but if rotate==0 we do nothing
          temp[0] = pdCode[components[i+1]][6];
          pdCode[components[i+1]][6] = pdCode[components[i+2]-1][6];
          pdCode[components[i+2]-1][6] = temp[0];
          int k = 1;
          for(int j = components[i+1]; j < (components[i+2]+components[i+1])/2; j++){
            swapTanglesAB(j, components[i+2]-k, newRow, pdCode, row, edge);
            k++;
          }
          for(int j = components[i+1]; j < components[i+2]; j++){
            for(int k = 0; k < 4; k++){
              temp[k] = pdCode[j][k];
            }for(int k = 0; k < 4; k++){
              pdCode[j][k] = temp[(k+2)%4];
            }
            for(int k = 0; k < 2*row; k++){
              if(edge[k][0]==j){
                edge[k][1] = (edge[k][1] + 2)%4;
              }
              if(edge[k][2]==j){
                edge[k][3] = (edge[k][3] + 2)%4;
              }
            }
          }
        }
        pdCode[components[i+1]-1][6] = 1;
        pdCode[components[i+2]-1][6] = 0;
        top_right = i+1;
      }
    }  
    else if((i==0 && pdCode[0][6]!=1)||previous_connection==-1){//then this is a mult component
      //Repeat first if *but* check for sum first, product second
      //a sum uses Clock 2 first tangle, and clock 1 last tangle
      
      int top_tangle = components[top_right];
      if(pdCode[components[top_right]][6]==1){
        top_tangle = components[top_right+1]-1;
      }
      int bottom_tangle = components[bottom_left+1]-1;
      
      getTangle2(row, top_tangle, 2, &tangle2i, &tangle2iClock, pdCode, edge);
      getTangle2(row, bottom_tangle, 1, &tangle2ii, &tangle2iiClock, pdCode, edge);

      int here = 0;
      for(int j = i; j < pieces ; j++){
        if(tangle2i < components[j+1] && tangle2ii < components[j+1]){
          here = j;
          break;
        }
      }
      
      if(here!=0 && ((tangle2i==tangle2ii)||(tangle2i==components[here] && tangle2ii==components[here+1]-1)||(tangle2i==components[here+1]-1 && tangle2ii==components[here]))){
      
        int rotate = tangle2iiClock; //0 is none, 1 is 90 CW, 2 is 180, 3 is 90 CCW
        if(here > i+1){//Then we need to swap these components
          int size = 0;
          int diff = 0;
          
          if(components[here+1] - components[here] >= components[i+2] - components[i+1]){
            size = components[i+2] - components[i+1];
            diff = components[here+1] - components[here] - size;
            insertZeroRows(components[i+2]-1, diff, newRow, pdCode);
            for(int j = components[i+1]; i < components[i+1]+size+diff; j++){
              swapTanglesAB(j, components[here]+j-components[i+1], newRow, pdCode, row, edge);
            }
            for(int j = components[here+1]-diff; j < newRow + diff; j++){
              for(int k = 0; k < 7; k++){
                pdCode[j][k] = pdCode[j + diff][k];
              }
            }
          } else {
            size = components[here+1] - components[here];
            diff = components[i+2] - components[i+1] - size;
            insertZeroRows(components[here+1]-1, diff, newRow, pdCode);
            for(int j = components[i+1]; j < components[i+1]+size+diff; j++){
              swapTanglesAB(j, components[here]+j-components[i+1], newRow, pdCode, row, edge);
            }
            
            for(int j = components[i+2]-diff; j < newRow + diff; j++){
              for(int k = 0; k < 7; k++){
                pdCode[j][k] = pdCode[j + diff][k];
                if((k < 4) && edge[pdCode[j][k]][0]==j+diff){
                  edge[pdCode[j][k]][0] -= diff;
                }
                if((k < 4) && edge[pdCode[j][k]][2]==j+diff){
                  edge[pdCode[j][k]][2] -= diff;
                }
              }
            }
          }
          
          int x = 1;
          for(int j = 1; j < newRow; j++){
            if(pdCode[j-1][6]==0){
              components[x] = j;
              x++;
            }
          }
        }
        
        if(rotate%2 == 1) {
          //0 is none, 1 is 90 CW, 2 is 180, 3 is 90 CCW
          //rotate = 1 => 0->3->2->1->0
          //rotate = 3 => 0->1->2->3->0
          for(int j = components[i+1]; j < components[i+2]; j++){
            rotateTangleFraction(pdCode, j);
            for(int k = 0; k < 4; k++){
              temp[k] = pdCode[j][k];
            }
            for(int k = 0; k < 4; k++){
              pdCode[j][k] = temp[(k+rotate)%4];
            }
            for(int k = 0; k < 2*row; k++){
              if(edge[k][0]==j){
                edge[k][1] = (edge[k][1] +rotate)%4;
              }
              if(edge[k][2]==j){
                edge[k][3] = (edge[k][3] +rotate)%4;
              }
            }
            pdCode[j][6] *= -1;
          }
          if((rotate == 1 && pdCode[components[i+1]][6] == 1)||(rotate == 3 && pdCode[components[i+1]][6]== -1)){
            temp[0] = pdCode[components[i+1]][6];
            pdCode[components[i+1]][6] = pdCode[components[i+2]-1][6];
            pdCode[components[i+2]-1][6] = temp[0];
            int k = 1;
            for(int j = components[i+1]; j < (components[i+2]+components[i+1])/2; j++){
              swapTanglesAB(j, components[i+2]-k, newRow, pdCode, row, edge);
              k++;
            }
            
          }
        }
        else if (rotate==2){
          temp[0] = pdCode[components[i+1]][6];
          pdCode[components[i+1]][6] = pdCode[components[i+2]-1][6];
          pdCode[components[i+2]-1][6] = temp[0];
          
          int k = 1;
          for(int j = components[i+1]; j < (components[i+2]+components[i+1])/2; j++){
            swapTanglesAB(j, components[i+2]-k, newRow, pdCode, row, edge);
            k++;
          }
          for(int j = components[i+1]; j < components[i+2]; j++){
            for(int k = 0; k < 4; k++){
              temp[k] = pdCode[j][k];
            }
            for(int k = 0; k < 4; k++){
              pdCode[j][k] = temp[(k+2)%4];
            }
            for(int k = 0; k < 2*row; k++){
              if(edge[k][0]==j){
                edge[k][1] = (edge[k][1] + 2)%4;
              }
              if(edge[k][2]==j){
                edge[k][3] = (edge[k][3] + 2)%4;
              }
            }
          }
        }
        pdCode[components[i+1]-1][6] = 1;
        pdCode[components[i+2]-1][6] = 0;
        top_right=i+1;
        printf("\n");
        
        continue;
      }
      //a product uses Clock 0 last tangle, and Clock 1 last tangle
      int left_tangle = components[bottom_left];
      int right_tangle = components[bottom_left+1]-1;
      if(pdCode[components[bottom_left]][6]==-1){
        left_tangle = right_tangle;
      }
      getTangle2(row, left_tangle, 0, &tangle2i, &tangle2iClock, pdCode, edge);
      getTangle2(row, right_tangle, 1, &tangle2ii, &tangle2iiClock, pdCode, edge);

      here = 0;
      for(int j = i; j < pieces; j++){
        if(tangle2i < components[j+1] && tangle2ii < components[j+1]){
          here = j;
          break;
        }
      }
      if(here!=0&&((tangle2i==tangle2ii)||(tangle2i==components[here] && tangle2ii==components[here+1]-1)||(tangle2ii==components[here] && tangle2i==components[here + 1] - 1))){
        //Then the tangle(s) belong to the same component so we can proceed
        int rotate = tangle2iiClock;//2 is none, 3 is 90 cw, 0 is 180, 1 is 90 ccw
        if(here > i+1){//Then we need to swap these components
          int size = 0;
          int diff = 0;
          if(components[here+1] - components[here] >= components[i+2] - components[i+1]){
            size = components[i+2] - components[i+1];
            diff = components[here+1]-components[here] - size;
            insertZeroRows(components[i+2]-1, diff, newRow, pdCode);
            //included zero rows to components[i+1] piece to make them equal sizes
            for(int j = components[i+1]; j < components[i+1] + size + diff; j++){
              swapTanglesAB(j, components[here] + j - components[i+1], newRow, pdCode,row, edge);
            }
            for(int j = components[here+1]-diff;j < newRow+diff; j++){
              for(int k = 0; k < 7; k++){
                pdCode[j][k] = pdCode[j+diff][k];
              }
            }
          } else {
            size = components[here+1] - components[here];
            diff = components[i+2] - components[i+1] - size;
            insertZeroRows(components[here+1]-1, diff, newRow, pdCode);
            //included zero rows to components[location] piece to make them equal sizes
          
            for(int j = components[i+1]; i < components[i+1] + size + diff; j++){
              swapTanglesAB(j, components[here] + j - components[i+1], newRow, pdCode,row, edge);
            }
            for(int j = components[i+2]-diff; j < newRow + diff; j++){
              for(int k = 0; k < 7; k++){
                pdCode[j][k] = pdCode[j+diff][k];
                if((k < 4) && edge[pdCode[j][k]][0]==j+diff){
                  edge[pdCode[j][k]][0] -= diff;
                }
                if((k < 4) && edge[pdCode[j][k]][2]==j+diff){
                  edge[pdCode[j][k]][2] -= diff;
                }
              }
            }
          } 
        //Components are swapped, one at location is now at i+1*******************
          //now update components array
          int x = 1;
          for(int j = 1; j < newRow; j++){
            if(pdCode[j-1][6]==0){
              components[x] = j;
              x++;
            }
          }
        }
        //Components are swapped and component array is updated
        //rotate component as needed to perform the operation
        if(rotate%2==1){//odd rotations need flipped fractions a/b -> -b/a
          for(int j = components[i+1]; j < components[i+1]-1; j++){
            rotateTangleFraction(pdCode, j);
            for(int k = 0; k < 4; k++){
              temp[k] = pdCode[j][k];
            }
            for(int k = 0; k < 4; k++){
              //rotate = 1 is 90 CCW (3->0->1->2->3); orig -> (orig + rotate)%4
              //rotate = 3 is 90 CW (0->3->2->1->0); orig -> (orig + rotate)%4
              pdCode[j][k] = temp[(k + 2 + rotate)%4];
            }
            for(int k = 0; k < 2*row; k++){
              if(edge[k][0]==j){
                edge[k][1] = (edge[k][1] + 2 + rotate)%4;
              }
              if(edge[k][2]==j){
                edge[k][3] = (edge[k][3] + 2 + rotate)%4;
              }
            }
            pdCode[j][6]*=-1;
          }
          if((rotate==1 && pdCode[components[i+1]][6] == -1) || (rotate == 3 && pdCode[components[i+1]][6] == 1)){//CCW makes first in sum, last in product
            temp[0] = pdCode[components[i+1]][6];
            pdCode[components[i+1]][6] = pdCode[components[i+2]-1][6];
            pdCode[components[i+2]-1][6] = temp[0];
            int k = 1;
            for(int j = components[i+1]; j < (components[i+2]+components[i+1])/2; j++){
              swapTanglesAB(j, components[i+2]-k, newRow, pdCode, row, edge);
              k++;
            }
          }
        }
        else if (rotate == 0){
          temp[0] = pdCode[components[i+1]][6];
            pdCode[components[i+1]][6] = pdCode[components[i+2]-1][6];
            pdCode[components[i+2]-1][6] = temp[0];
          int k = 1;
          for(int j = components[i+1]; j < (components[i+2]+components[i+1])/2; j++){
            swapTanglesAB(j, components[i+2]-k, newRow, pdCode, row, edge);
            k++;
          }
          for(int j = components[i+1]; j < components[i+2]; j++){
            for(int k = 0; k < 4; k++){
              temp[k] = pdCode[j][k];
            }
            for(int k = 0; k < 4; k++){
              pdCode[j][k] = temp[(k + 2)%4];
            }
            for(int k = 0; k < 2*row; k++){
              if(edge[k][0]==j){
                edge[k][1] = (edge[k][1] + 2)%4;
              }
              if(edge[k][2]==j){
                edge[k][3] = (edge[k][3] + 2)%4;
              }
            }
          }
        }
        pdCode[components[i+1]-1][6] = -1;
        pdCode[components[i+2] - 1][6] = 0;
        bottom_left = i+1;
      }
      
    }
    else{//This is a single rational tangle in combination with montesinos
      //Repeat either if, but sums and products use only single tangle.
      //might need to check how single component combines with the previous one
      printf("Possibly non-algebraic, N(");
      for(int i = 0; i < newRow-1; i ++){
        char op = '?';
        if(pdCode[i][6]==1){
          op = '+';
        } else if (pdCode[i][6]==-1){
          op = '*';
        }
        printf("%d/%d %c", pdCode[i][4], pdCode[i][5], op);
      }
      printf("%d/%d),,,,", pdCode[newRow-1][4], pdCode[newRow-1][5]);
      exit(0);
    }
  }
}

/*
  Should analyze the non-combinable tangles and produce a continued
  fraction decomposition in the form of a sum of continued
  fractions N(a_1/b_1 + ... a_n/b_n).
*/
void algTangle(int row, int newRow, int pdCode[][7], int edge[][4]){ 
  int tangle2i, tangle2iClock, tangle2ii, tangle2iiClock;
  sort(row, newRow, pdCode, edge);
  int temp[4];
  int components[newRow];
  for(int i = 0; i < newRow; i++){
    components[i]=0;
  }
  
  //Take pdCode and add new column for operations
  //Then tangle at i attaches to i+1 using operation pdCodePlus[i][7]
  //-1 = *, 0 = ?, 1 = +

  int rotated=0;
  int start = 0;
  int Tangle = 0;
  int adding = 0;
  int stuck = 0;
  int k = 1;
  while (Tangle < newRow-1){
    
    //look for horizontal sums along right side
    while(adding == 0){
      getTangle2(row, Tangle, 1, &tangle2i, &tangle2iClock, pdCode, edge);
      getTangle2(row, Tangle, 2, &tangle2ii, &tangle2iiClock, pdCode, edge);
      if(tangle2i == tangle2ii && tangle2i > Tangle && pdCode[tangle2i][4] != 0){
        stuck = 0;
        //storing found tangle
        for(int i = 0; i < 4; i++){
          temp[i] = pdCode[tangle2i][i];
        }
        /*
          In a straight horizontal product, tangle2iClock = 0 so 
            (i + tangle2iClock)%4 = (i + 0)%4 = i
            reflecting no need for pdCode shift.
          In horizontal product with ccw rotation, tangle2iClock = 3 so
            (i + tangle2iClock)%4 = (i + 3)%4 
            reflecting a right shift of pdCode.
          In horizontal product with 180 rotation, tangle2iClock = 2 so
            (i + tangle2iClock)%4 = (i + 2)%4
            reflecting a two space shift of pdCode.
          In horizontal product with cw rotation, tangle2iClock = 1 so
            (i + tangle2iClock)%4 = (i + 1)%4
            reflecting a left shift of pdCode.
        */
        for(int i = 0; i < 4; i++){
          pdCode[tangle2i][i] = temp[(i + tangle2iClock) % 4];
        }
        int clockAdjust = 0;
        if(tangle2iClock%2 == 0){  
          clockAdjust = tangle2iClock;
        } 
        else {
          clockAdjust = (tangle2iClock + 2)%4;
          rotateTangleFraction(pdCode, tangle2i);
        }
        for(int i = 0; i < 2*row; i++){
          if(edge[i][0] == tangle2i){
            edge[i][1] = (edge[i][1] + clockAdjust)%4;
          }
          if(edge[i][2] == tangle2i){
            edge[i][3] = (edge[i][3] + clockAdjust)%4;
          }
        }
        if(tangle2i > Tangle+1){
          swapTanglesAB(Tangle+1, tangle2i, newRow, pdCode, row, edge);
        }
        pdCode[Tangle][6] = 1;
        Tangle++;
      }
      else {
        adding = 1;
        stuck++;
      }
    }
    //look for vertical sums along bottom side
    while(adding == 1){      
      getTangle2(row, Tangle, 0, &tangle2i, &tangle2iClock, pdCode, edge);
      getTangle2(row, Tangle, 1, &tangle2ii, &tangle2iiClock, pdCode, edge);
      if(tangle2i == tangle2ii && tangle2i > Tangle && pdCode[tangle2i][4] != 0){
        stuck = 0;
        for(int i = 0; i < 4; i++){
          temp[i] = pdCode[tangle2i][i];
        }
        /*
          In a straight vertical product, tangle2iiClock = 2 so 
            (i + tangle2iiClock + 2)%4 = (i + 4)%4 = i
            reflecting no need for pdCode shift.
          In vertical product with ccw rotation, tangle2iiClock = 1 so
            (i + tangle2iiClock + 2)%4 = (i + 3)%4 
            reflecting a right shift of pdCode.
          In vertical product with 180 rotation, tangle2iiClock = 0 so
            (i + tangle2iiClock + 2)%4 = (i + 2)%4
            reflecting a two space shift of pdCode.
          In vertical product with cw rotation, tangle2iiClock = 3 so
            (i + tangle2iiClock + 2)%4 = (i + 5)%4 = (i + 1)%4
            reflecting a left shift of pdCode.
        */
        for(int i = 0; i < 4; i++){
          pdCode[tangle2i][i] = temp[(i + tangle2iiClock + 2) % 4];
        }
        int clockAdjust = 0;
        if(tangle2iClock%2 == 0){  
          clockAdjust = tangle2iiClock;
          rotateTangleFraction(pdCode, tangle2i);
        } 
        else {
          clockAdjust = (tangle2iiClock + 2)%4;
        }
        for(int i = 0; i < 2*row; i++){
          if(edge[i][0] == tangle2i){
            edge[i][1] = (edge[i][1] + clockAdjust)%4;
          }
          if(edge[i][2] == tangle2i){
            edge[i][3] = (edge[i][3] + clockAdjust)%4;
          }
        }
        if(tangle2i > Tangle+1){
          swapTanglesAB(Tangle+1, tangle2i, newRow, pdCode, row, edge);
        }
        pdCode[Tangle][6] = -1;
        Tangle++;
        
      }
      else {
        adding = 0;
        stuck++;
        //Tangle++;
      } 
    }
    if(stuck > 1){
      if(start == Tangle && rotated==0){
        //try rotating the start piece by 180
        for(int i = 0; i < 4; i++){
          temp[i] = pdCode[start][i];
          if(edge[temp[i]][0]==Tangle){
          edge[temp[i]][1] = (edge[temp[i]][1] + 2)%4;
          } else if (edge[temp[i]][2]==Tangle){
          edge[temp[i]][3] = (edge[temp[i]][3] + 2)%4;
          }
        }
        for(int i = 0; i < 4; i++){
          pdCode[start][i] = temp[(i+2)%4];
        }
        stuck=0;
        rotated = 1;
      } else {
        Tangle++;
        components[k] = Tangle;
        k++;
        start = Tangle;
        rotated = 0;
      }
    }
  }
  
  if(components[k-1]!= newRow && pdCode[newRow - 1][6] == 0){
    components[k] = newRow;
    k++;
  }
  
  if(/*stuck > 1 && */k > 2){
  
    k--;
    
    
    //Orient and order of components for algebraic
    // if K is even, first piece should be oppositely oriented
      //i.e. expected mont operation is * or (-1)
    // if k is odd, first piece should be normally oriented
      //i.e. expected mont operation is + or (1)
    //up to here, pdCode, components, and edge should be accurate, now to verify operations
    
    orientAlgebraic(k, newRow, row, components, pdCode, edge);
    
    //Make sure right most component is summed to the rest, if not, rotate by 90
    if(pdCode[components[k-1]-1][6]== -1){//Then rotate everything 90 CCW Clocks 0->1->2->3->0, op 1 <=> -1
      //When turning product components into summed components, tangle order remains
      //When turning summed components into product components, tangle order *reverses*

      for(int i = 0; i < k; i++){
        for(int j = components[i]; j < components[i+1]; j++){
          rotateTangleFraction(pdCode, j);
          pdCode[j][6] *= -1;
          for(int l = 0; l < 4; l++){
            temp[l] = pdCode[j][l];
          }
          for(int l = 0; l < 4; l++){
            pdCode[j][l] = temp[(l+3)%4];
          }
          for(int l = 0; l < 2*row; l++){
            if(edge[l][0]==j){
              edge[l][1] = (edge[l][1] + 3)%4;
            }
            if(edge[l][2]==j){
              edge[l][3] = (edge[l][3] + 3)%4;
            }
          }
        }
        if(pdCode[components[i]][6]== -1){//Then this component became a product comp, so was a summed comp and tangle order reverses
          temp[0] = pdCode[components[i]][6];
          pdCode[components[i]][6] = pdCode[components[i+1]-1][6];
          pdCode[components[i+1]-1][6] = temp[0];
          int fromEnd = 1;
          for(int j = components[i]; j < (components[i+1]+components[i])/2; j++){
            swapTanglesAB(j, components[i+1]-fromEnd, newRow, pdCode, row, edge);
            fromEnd++;
          }
        }
      }
    }

    for(int i = 0; i < newRow; i++){
      if(pdCode[i][5] < 0){
        pdCode[i][4] *= -1;
        pdCode[i][5] *= -1;
      }
    }
     

    printf("N(");    
    for(int i = 0; i < k; i++){
      int mont=isMontesino(components[i], components[i+1], pdCode, 1);
      handleMontesinosComponent(components[i], components[i+1],pdCode, mont);
      if(i < k-1 && pdCode[components[i+1]-1][6]==1){  
        printf("+");
      }
      if(i < k-1 && pdCode[components[i+1]-1][6]==-1){
        printf("*");
      }
    }

    printf("), N(");
    int remainder[k];
    int sign[k];
    for(int i = 0; i < k; i++){
      remainder[i] = 0;
      sign[i] = 0;
    }
      
    for(int i = 0; i < k; i++){
      
      sign[i] = majoritySign(components[i], components[i+1], pdCode, &remainder[i], k-i);
      makeCanonical(components[i], components[i+1], pdCode, sign[i], &remainder[i], k - i);
      //Might need to split getFrac from the rest so I can pull common negative to the front
      getFraction(components[i], components[i+1], newRow, pdCode, &remainder[i], sign[i], k-i);
      if(i < k-1){
        if(pdCode[components[i+1] - 1][6] == 1){
          printf(" + ");
        } else if (pdCode[components[i+1] - 1][6] == -1){
          printf(" * ");
        }
      }
    }
    printf("),[");
    for(int i = 0; i < k; i++){
      getConwayMontesinos(sign[i], remainder[i], components[i], components[i+1], newRow, pdCode);
      if(i < k-1){
        printf("|");
      }
    }
    printf("],,,");
\
  } else {
    //standard Montesinos case
    handleMontesinos(newRow, pdCode);
  }
}

// Given Tangle A, finds a second tangle which can be combined with it along the side specified to create
// a horizontal or vertical tangle
int makeHorVtangle(int TangleA, int Clock1, int Clock2, int rows, int pdCode[][7], int edge[][4]){
  int TangleB = 0;
  int rotate = 0;
  int arc1, arc2, tangle2i, tangle2ii, tangle2iClock, tangle2iiClock;

  if(pdCode[TangleA][4] != 0 && ( abs(pdCode[TangleA][4]) == 1 || pdCode[TangleA][5] == 1 ) ){ 

    TangleB = canCombine(TangleA, Clock1, Clock2, rows, pdCode, edge);
      
    if(TangleB > -1 && TangleB < TangleA){
      arc1 = pdCode[TangleA][Clock1];
      arc2 = pdCode[TangleA][Clock2];
    
      //Get clock positions of arcs in TangleB to rotate as needed
      getTangle2(rows, TangleA, Clock1, &tangle2i, &tangle2iClock, pdCode, edge);
      getTangle2(rows, TangleA, Clock2, &tangle2ii, &tangle2iiClock, pdCode, edge);
      
      if(tangle2iClock != (Clock1 + 3) % 4){
        //Need to rotate TangleB
        rotate = (Clock1 + 3 - tangle2iClock) % 4;
        rotateTangle90CCW_Ntimes(rotate, TangleB, rows, pdCode, edge);
      }
    
    //printf("Combinging TangleA = %d to TangleB = %d after rotating TangleB %d along edge %d %d\n", TangleA, TangleB, rotate, Clock1, Clock2);
      if((Clock1 % 2) == 0 && abs(pdCode[TangleB][4]) == 1 && abs(pdCode[TangleA][4]) == 1){
        // Clock1 even => vertical sum, need unit numerators
        
        if(pdCode[TangleB][4] < 0){
          pdCode[TangleB][4] *= -1;
          pdCode[TangleB][5] *= -1;
        }
        if(pdCode[TangleA][4] < 0){
          pdCode[TangleA][4] *= -1;
          pdCode[TangleA][5] *= -1;
        }
        pdCode[TangleB][5] += pdCode[TangleA][5];

        if(pdCode[TangleB][5] < 0){
          pdCode[TangleB][4] *= -1;
          pdCode[TangleB][5] *= -1;
        }      
        
        combineTanglesPDandEdge(TangleA, TangleB, Clock1, Clock2, rows, pdCode, edge);
        
        return 1;
      } else if((Clock1 % 2) == 1 && pdCode[TangleA][5] == 1 && pdCode[TangleB][5] == 1){
        //Clock1 odd => horizontal sum, need denom of 1
        pdCode[TangleB][4] += pdCode[TangleA][4];
        
        combineTanglesPDandEdge(TangleA, TangleB, Clock1, Clock2, rows, pdCode, edge);
        
        return 1;
      } else {
        return 0;
      }
    }
    else {
      return 0;
    }
  } else {
    return 0;
  }
}

int makeRationalSimple(int TangleA, int Clock1, int Clock2, int newRow, int rows, int pdCode[][7], int edge[][4]){
  if(pdCode[TangleA][4] == 0){
    return 0;
  }
  int TangleB = 0;
  int rotate = 0;
  int tangle2i, tangle2ii, tangle2iClock, tangle2iiClock;
  int arc1 = pdCode[TangleA][Clock1];
  int arc2 = pdCode[TangleA][Clock2];

  TangleB = canCombine(TangleA, Clock1, Clock2, rows, pdCode, edge);
  if(TangleB > -1){
    
    //find where arc1 and arc2 connect to TangleB to rotate if needed
    getTangle2(rows, TangleA, Clock1, &tangle2i, &tangle2iClock, pdCode, edge);
    getTangle2(rows, TangleA, Clock2, &tangle2ii, &tangle2iiClock, pdCode, edge);
    
    if(tangle2iClock != (Clock1 + 3) % 4){
      rotate = (Clock1 + 3 - tangle2iClock) % 4;
      rotateTangle90CCW_Ntimes(rotate, TangleB, rows, pdCode, edge);
    }
    
    //Need to know if this is a vertical or horizontal combination
    if( Clock1 % 2 == 1 && ( pdCode[TangleB][5] == 1 || pdCode[TangleA][5] == 1 ) ){
      
      //This case is a horizontal combination with horizontal tangle
      int horizontal = TangleB;
      int rational = TangleA;
      if(pdCode[TangleB][5] != 1){
        rational = TangleB;
        horizontal = TangleA;
      }
      //If TangleA has fraction a/b and TangleB has fraction n/1
      //Write a = a'b + r, then a/b + n/1 = r/b + a' + n 
      // = (r + b(a' + n)) / b  is the resulting fraction
      // a' is the same as floor(a/b) and r is the same as a modulo b 
      
      int a = pdCode[rational][4];
      int b = pdCode[rational][5];
      int n = pdCode[horizontal][4];

      pdCode[TangleB][4] = a + b*n;
      pdCode[TangleB][5] = b;
      
      combineTanglesPDandEdge(TangleA, TangleB, Clock1, Clock2, rows, pdCode, edge);
      
      return 1;
    } else if( Clock1 % 2 == 0 && ( abs(pdCode[TangleB][4]) == 1 || abs(pdCode[TangleA][4]) == 1 ) ) {
      //This case is vertical combination with vertical tangle
      int rational = TangleA;
      int vertical = TangleB;

      if( abs(pdCode[TangleB][4]) != 1){
       rational = TangleB;
       vertical  = TangleA; 
      }
      //Push negative to denominator if needed
      if(pdCode[vertical][4] < 0){
        pdCode[vertical][4] *= -1;
        pdCode[vertical][5] *= -1;
      }
      
      //If TangleA has fraction a/b and TangleB has fraction 1/n
      //Resulting  fraction is 1/ (n + 1/(a/b)) = 1/(n + b/a) 
      // or a/( a*n + b )
      int a = pdCode[rational][4];
      int b = pdCode[rational][5];
      int n = pdCode[vertical][5];

      pdCode[TangleB][4] = a;
      pdCode[TangleB][5] = a*n + b;

      if(pdCode[TangleB][5] < 0){
        pdCode[TangleB][4] *= -1;
        pdCode[TangleB][5] *= -1;
      }
      combineTanglesPDandEdge(TangleA, TangleB, Clock1, Clock2, rows, pdCode, edge);
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

void pdToConway(int row, int pdCode[][7]){

  int edge[2 * row + 2][4];
  int i, j, newRow, Tangle, signCrossing[row], crossing2, crossing2Clock,
      writhe, temp;

  /* Create edge matrix where each ROW corresponds to an Arc */
  createEdge(row, pdCode, edge);
  
  writhe = compute_writhe(row, pdCode, edge);
  
  printf("%d,", writhe);

  /*******************************************************************/
  /**** combine 1/1 tangles into (n/1) and (1/n) tangles  **************/
  /*******************************************************************/

  /* 
   Combine 1/1 tangles into horizontal (n/1) integer tangles: Old + New -> New
   Cycling through pdCode from bottom to top, ignores first row.
   */

  /*
  First we go through to add another tangle to the right of the current, resulting 
  tangle row will have 0, 1 arc from Tangle then 2, 3 arc from new tangle. Unless 
  rotating is needed. Then 0, 1 arc from Tangle then 2, 3 arc from new tangle's 0, 1.
  */

  int canADD,arc1, arc2, tangle2i, tangle2ii, tangle2iClock, tangle2iiClock;
  int added = 1;

  //MAKING A VERTIC SUM INTO HORIZONTAL FRAC. ROWS 1 AND 3 MAKE -1/2 NOT -2/1
  while (added > 0){
    added = 0;
    Tangle = row - 1;
    //make simple by adding basic
    while(Tangle > 0){
      for(int i = 0; i < 4; i++){
        added +=makeHorVtangle(Tangle, i, (i + 1) % 4, row, pdCode, edge);
      }
      Tangle--;
    }
  }
  
  

  newRow = removeTangles(row, pdCode, edge, row);
  added = 1;

  while(added > 0){
    //make rationals by adding simple
    Tangle = newRow - 1;
    added = 0;
    while(Tangle > 0){      
      for(int i = 0; i < 4; i++){
        added += makeRationalSimple(Tangle, i, (i + 1) % 4, newRow, row, pdCode, edge);
      }
      Tangle--;
    }
  }
  
  newRow = removeTangles(row, pdCode, edge, newRow);
  
  if(newRow == 2){
    addRationalTangles(row, pdCode, edge);
    newRow = removeTangles(row, pdCode, edge, newRow);
  }


  
    
    // while (Tangle > 0) {
    //   if (pdCode[Tangle][4] != 0){
    //     canADD = addTangles(row, Tangle, 1, 2, pdCode, edge);
    //     if (canADD > 0){
    //       added = 1;
    //       printf("Right sum to %d\n", Tangle);
    //       display(newRow, 7, pdCode);
    //     }
    //   }
    //   Tangle = Tangle - 1;
    // }
    // /* 
    // combine 1/1 tangles into horizontal (n/1) integer tangles:
    //  New + Old -> New
    // */
    // Tangle = row - 1;
    // while (Tangle > 0) {
    //   if (pdCode[Tangle][4] != 0){
    //     canADD = addTangles(row, Tangle, 3, 0, pdCode, edge);
    //     if (canADD > 0){
    //       added = 1;
    //       printf("Left sum to %d\n", Tangle);
    //       display(newRow, 7, pdCode);
    //     }
    //   }
    //   Tangle = Tangle - 1;
    // }

    // /* combine 1/1 tangles into vertical (1/n) tangles: Old * New -> New  */
    // Tangle = row - 1;
    // while (Tangle > 0) {
    //   if (pdCode[Tangle][4] != 0){
    //     canADD = addTangles(row, Tangle, 1, 0, pdCode, edge);
    //     if (canADD > 0){
    //       added = 1;
    //       printf("Bottom sum to %d\n", Tangle);
    //       display(newRow, 7, pdCode);
    //     }
    //   }
    //   Tangle = Tangle - 1;
    // }
    
    // /* combine 1/1 tangles into vertical (1/n) tangles: New * Old -> New  */
    // Tangle = row - 1;
    // while (Tangle > 0) {
    //   if (pdCode[Tangle][4] != 0){
    //     canADD = addTangles(row, Tangle, 3, 2, pdCode, edge);
    //     if (canADD > 0){
    //       added = 1;
    //       printf("Top sum to %d\n", Tangle);
    //       display(newRow, 7, pdCode);
    //     }
    //   }
    //   Tangle = Tangle - 1;
    // }
    /* remove rows in pdCode where tangle has been combined with another tangle
    This statement feels unnecessary newRow = row;, so far its not causing issues to remove.
    newRow here adjusts pdCode and edges to move zero rows to the bottom of pdCode
    then relabel crossing numbers in edge to match the new pdCode.
    newRow then holds the location of the last non-zero row.
    */
    // newRow = removeTangles(row, pdCode, edge, row);
    // /*******************************************************************/
    // /****** combine (n) and (1/n) tangles into rational tangles  *******/
    // /*******************************************************************/

    // // perform horizontal sums
    // Tangle = newRow - 1;
    // while (Tangle > 0) {
    //   if (pdCode[Tangle][4] != 0){
    //     canADD = addTangles(row, Tangle, 3, 0, pdCode, edge);
    //     if (canADD > 0){
    //       added = 1;
    //     }
    //   }
    //   Tangle = Tangle - 1;
    // }

    // Tangle = newRow - 1;
    // while (Tangle > 0) {
    //   if (pdCode[Tangle][4] != 0){
    //     canADD = addTangles(row, Tangle, 1, 2, pdCode, edge);
    //     if (canADD > 0){
    //       added = 1;
    //     }
    //   }
    //   Tangle = Tangle - 1;
    // }

    // // done with horizontal sums
    // //***************************************
    // // perform vertical sums (i.e., products)
    // Tangle = newRow - 1;
    // while (Tangle > 0) {
    //   if (pdCode[Tangle][4] != 0){
    //     canADD = addTangles(row, Tangle, 1, 0, pdCode, edge);
    //     if (canADD > 0){
    //       added = 1;
    //     }
    //   }
    //   Tangle = Tangle - 1;
    // }
    // Tangle = newRow - 1;
    // while (Tangle > 0) {
    //   if (pdCode[Tangle][4] != 0){
    //     canADD = addTangles(row, Tangle, 3, 2, pdCode, edge);
    //     if (canADD > 0){
    //       added = 1;
    //     }
    //   }
    //   Tangle = Tangle - 1;
    // }

    

    // // done with vertical sums
    // /*****************************************************************/
    // /* 
    //   Adjust edge and pdCode to move zero rows to the bottom and relabel
    //   crossings in edge to accomodate. Then save the location of the last non-zero
    //   row in pdCode.
    // */
    
    // newRow = removeTangles(row, pdCode, edge, newRow);
    // /*
    //   Entering this conditional statement occurs if the pdCode passes through
    //   the above attempts to add tangles, but can not make anymore integer tangles
    //   leaving behind 3 or more tangles which cannot be added to form a single
    //   rational tangle.
    // */

    if (added == 0 && newRow > 2){
      
      algTangle(row, newRow, pdCode, edge);
      exit(0);
    }
  
  /* 
    If newRow = 2, there are only 2 rational tangles which should combine to another rational knot. 
   */
  // if (newRow == 2) {
  //   if (pdCode[0][0] == pdCode[1][0]) {
  //     // If both tangles share a northwest arc, must rotate on of them 90 degrees.
  //     if (!(pdCode[0][1] == pdCode[1][3] && pdCode[0][2] == pdCode[1][2] &&
  //         pdCode[0][3] == pdCode[1][1])){
  //           //printf("Last two tangles are mismatched\n");
  //         }
  //     // Reflects the rotation in the fraction representation of the second tangle.
  //     temp = pdCode[1][4];
  //     pdCode[1][4] = pdCode[1][5];
  //     pdCode[1][5] = -temp;
  //     /**********************************************
  //       This doesn't rotate the actual pdCode though?
  //      **********************************************/
  //   } else if (pdCode[0][0] == pdCode[1][1]) {
  //     if (pdCode[0][1] != pdCode[1][0] || pdCode[0][2] != pdCode[1][3] ||
  //         pdCode[0][3] != pdCode[1][2])
  //       DEBUG_PRINTF("MISTAKE1\n");
  //   } else if (pdCode[0][0] == pdCode[1][2]) {
  //     if (pdCode[0][1] != pdCode[1][1] || pdCode[0][2] != pdCode[1][0] ||
  //         pdCode[0][3] != pdCode[1][3])
  //       DEBUG_PRINTF("MISTAKE2\n");
  //     // rotate tangle1
  //     temp = pdCode[1][4];
  //     pdCode[1][4] = pdCode[1][5];
  //     pdCode[1][5] = -temp;

  //   } else if (pdCode[0][0] == pdCode[1][3]) {
  //     if (pdCode[0][1] != pdCode[1][2] || pdCode[0][2] != pdCode[1][1] ||
  //         pdCode[0][3] != pdCode[1][0])

  //       DEBUG_PRINTF("MISTAKE3\n");
  //   }
  //   newRow = addRationalTangles(row, pdCode, edge);
  // }

  // // newROW = 1 iff tangle is rational
  // if (newRow == 1) {
  //   if (pdCode[0][0] == pdCode[0][1] && pdCode[0][2] == pdCode[0][3]) {
  //     // rotate tangle0
  //     rotateTangleFraction(pdCode, 0);
  //   } else if (pdCode[0][0] != pdCode[0][3] || pdCode[0][1] != pdCode[0][2])
  //     DEBUG_PRINTF("mistake\n");
  // }
  

  /*******************************************************************/
  /****** if not rational, check for algebraic sums/products   *******/
  /*******************************************************************/

  /* if more than 2 rational tangles, create matrix with algebraic sums/products
   */
  // add this to:   addTangles(Tangle, ?, ?, pdCode,edge, AlgCode);
  // add new matrix AlgCode

  char add = '&';  // horizontal addition
  char mult = '*'; // vertical addition
  char divide = '/';
  char string[12]; // create an empty string to store number
  sprintf(string, "%d", pdCode[0][4]);
  strncat(string, &divide, 1);
  
  /* 
  Ugly looking manipulations
  Adjust fraction by rotating and cancelling negatives
  Creates a list of all equivalent (mirror and non-mirror) 
    fraction denominators, then outputs the smallest such in both cases.
  Output formatted to fill a .csv file via test.sh
  */
  //printf("N(%d/%d),", pdCode[0][4], pdCode[0][5]);

  // adjusting notation to match two-bridges a bad way
  
  int num;
  int den;
  
  if(pdCode[0][1]==0 || pdCode[0][3]==0){
    rotateTangleFraction(pdCode, 0);
  }
  if(pdCode[0][5] < 0){
    pdCode[0][4] *= -1;
    pdCode[0][5] *= -1;
  }
  int sign = 1;
  char pm = ' ';
  num = pdCode[0][4];
  den = pdCode[0][5];
  if(num < 0){
    sign  = -1;
    pm = '-';
    num *= -1;
  }
  printf("%d/%d,%c[", sign*num, den, pm);
  getConway(num, den);
  printf("],");
  sign *= -1;
  if(pm == '-'){
    pm = ' ';
  } else {
    pm = '-';
  }
  printf("%d/%d,%c[",sign*num, den, pm);
  getConway(num, den);
  printf("],");
  /* Routine for moving fraction to minimal in the context of knots and links
  int denList[abs(num)];
  i=0;
  int mirrorList[abs(num)];
  j = 0;
  if(abs(num) < abs(den)){
    int temp = num;
    num = -den;
    den = temp;
  }
  
  // Cancel double negative -p/-q = p/q
  if (num < 0 && den < 0){
    num = -num;
    den = -den;
  }

  // Write -p/q = p/-q
  if (num < 0 && den > 0){
    num = -num;
    den = -den;
  }
  
  printf("N(%d/%d),", num, den);
  
  denList[i] = aModB(den, num);
  ++i;

  int x = 1;
  while(x < num){
    if (aModB(x*den, num)==1 && x!=aModB(den, num)){
      denList[i] = x;
      ++i;
    }
    ++x;
  }

  mirrorList[j] = aModB(-den, num);
  ++j;

  x=1;
  while (x < num){
    if (aModB(x*den, num)==num-1 && x!=aModB(-den, num)){
      mirrorList[j] = x;
      ++j;
    }
    ++x;
  }

  x=0;
  temp=1;
  while (temp < i){
    if(denList[x] >= denList[temp]){
      x = temp;
    }
    ++temp;
  }
  
  printf("N(%d/%d),", num, denList[x]);
  
  //int conway[10] = {0,0,0,0,0,0,0,0,0,0};
  int a,b,q,k=0;
  a = num;
  b = denList[x];
  printf("[");
  getConway(a, b);
  printf("],");
  
  x = 0;
  temp = 1;
  while (temp < i){
    if (mirrorList[x] >= mirrorList[temp]){
      x = temp;
    }
    ++temp;
  }
  
  printf("N(%d/%d),", num, mirrorList[x]);
  a,b,q,k=0;

  a = num;
  b = mirrorList[x];
  
  printf("[");
  getConway(a,b);
  printf("],");
  */
}