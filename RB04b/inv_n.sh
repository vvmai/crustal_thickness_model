#!/bin/csh -f
#
# inversion of Rutter & Brodie (2004b) data on diffusion creep

#set up execution
set EXEC = labmc

set data = rutter04b_mini.dat
set outdir = .
set out = rutter04b_mini.out

echo $data

#execute code
#set conjugate = "-C10"
set conjugate = ""

$EXEC -D$data $conjugate -d1 -M10000/100/75/10 -V -R1 \
      -Pa-0/0/1/4/0/600/0/0  > $outdir/$out &

sleep 2

end