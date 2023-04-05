#!/bin/bash

#Test for -2/1 tangle
#printf "expecting tangle fraction -2/1 conway [-2]\n"
#./pdToConwayTangles "[[5,2,6,1],[2,5,3,4]]"

#Test below is for -5/3 tangle
#printf "expecting tangle fraction -5/3 conway -[2 1 1]\n"
#./pdToConwayTangles "[[1,5,2,4],[3,7,4,8],[5,10,6,9],[8,2,9,3]]"

#Test below is for -5/3 tangle
#printf "expecting tangle fraction -3/5 conway -[2 1 1 0]\n"
#./pdToConwayTangles "[[1,10,2,9],[3,7,4,8],[5,2,6,3],[8,4,9,5]]"

#Test below is for the 5/8 tangle
printf "expecting tangle fraction 5/8 conway [2 1 1 1 0]\n"
./pdToConwayTangles "[[6,2,7,1],[9,4,8,3],[11,8,10,7],[5,10,4,9],[2,12,3,11]]"
