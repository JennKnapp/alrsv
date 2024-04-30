#!/bin/bash

for  lab in "$@"; do
	alcov find_lineages --min_depth=5 --unique=False --csv=True --save_img=True "$lab".txt lineages.txt
done

for lab in "$@"; do
        alcov find_mutants --min_depth=5 --save_img=True "$lab".txt BA.4_BA.5_BA.5.2.1_BE.1.txt
done

