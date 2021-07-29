#!/bin/csh -f
#
# inversion of Fukuda (2018) data on wet composite creep

#set up execution
set EXEC = labmc

set data = fukuda18.dat
set out = fukuda18_f.out

echo $data

#execute code
set conjugate = ""

$EXEC -D$data $conjugate -d1 -d2 -M1000000/100/75/10 -V -R1 -B-0/0 \
      -Pf-0/0/1/6/1/6/1/5/0/600/0/0 > ./$out &

sleep 2