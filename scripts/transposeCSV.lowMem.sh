#!/bin/bash

# Script to transpose a csv file. Uses less memory, slower

numc=$(($(head -n 1 "$1" | grep -o , | wc -l)+1))
for ((i=1; i<="$numc"; i++))
do cut -d,  -f"$i" "$1" | paste -s -d,
done
