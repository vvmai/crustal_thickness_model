#!/bin/csh -f
#
# inversion of Luan & Paterson (1992) data on wet dislocation creep

#set up execution
set EXEC = labmc

set data = LP92.dat
set out = LP92_n.out

echo $data

#execute code
#set conjugate = "-C10"
set conjugate = ""

$EXEC -D$data $conjugate -d1 -d2 -M1000000/100/75/10 -V -R1 -B-5/5 \
      -Pd-0/0/1/6/1/5/0/600/0/0  > ./$out &

sleep 2

end