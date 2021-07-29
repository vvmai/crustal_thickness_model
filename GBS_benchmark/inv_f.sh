#!/bin/csh -f
#
# GBS inversion of Fukuda (2018) 

#set up execution
set EXEC = labmc

set data = syn_d.dat
set out = syn_d.out

echo $data

#execute code
set conjugate = ""

$EXEC -D$data $conjugate -d1 -d2 -M100000/100/75/10 -V -R1 -B-0/0 \
      -Pf-0/0/0/6/0/6/0/5/0/600/0/0 > ./$out &

sleep 2