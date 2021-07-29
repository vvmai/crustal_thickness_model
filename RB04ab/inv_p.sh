#!/bin/csh -f
#
# inversion of Rutter & Brodie (2004 a & b) data on diffusion and dislocation creep

#set up execution
set EXEC = labmc

set data = rutter04ab.dat
set out = rutter04ab_X.out

echo $data

set conjugate = ""

$EXEC -D$data $conjugate -d1 -M1000000/100/75/10 -V -R1 -B-5/5 \
      -Pa-0/0/1/6/0/600/0/0 -Pb-0/0/1/6/0/600/0/0 > ./$out &
      
sleep 2

end