#!/bin/csh -f
#
# inversion of Luan & Paterson (1992) data on wet dislocation creep

#set up execution
set EXEC = labmc

set data = ~/lab/LP92/LP92.dat
set out = LP92_pXC.out
set OUT_DIR = ~/lab/LP92

echo $data

#execute code
#set conjugate = "-C10"
set conjugate = "-C10"

$EXEC -D$data $conjugate -d1 -d2 -M1000000/100/75/10 -V -R1 -B-5/5 \
      -Pc-0/0/1/6/1/5/0/600/0/0 -Pd-0/0/1/6/1/5/0/600/0/0 > $OUT_DIR/$out &

sleep 2

end