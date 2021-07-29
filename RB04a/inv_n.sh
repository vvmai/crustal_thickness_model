#!/bin/csh -f
#
# inversion of Rutter & Brodie (2004a) data on dislocation creep

#set up execution
set EXEC = labmc

set data = rutter04a_1.dat
set outdir = .
set out = rutter04a_1.out

echo $data

#execute code
#set conjugate = "-C10"
set conjugate = ""

$EXEC -D$data $conjugate -d1 -M100000/100/75/10 -V -R1 \
      -Pb-0/0/1/6/0/600/0/0  > $outdir/$out &

sleep 2

end