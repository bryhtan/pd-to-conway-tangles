#!/bin/bash
file="pdCodes.txt"

echo "PD to Conway Tangles Results;" > tangle_out.csv

while IFS=$'\t' read -r twist frac pd
do
    [ -z "$pd" ] && continue;
    printf "%s\n" "$frac"
    printf "%s,%s," "$twist" "$frac">>tangle_out.csv
    ./pdToConwayTangles "$pd" >>tangle_out.csv
    echo ";" >>tangle_out.csv
done <$file
