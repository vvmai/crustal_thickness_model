#!/bin/csh -f
#
# inversion of Gleason & Tullis (1995) data on dry dislocation creep

#set up execution
set EXEC = labmc

set data = GT95_d.dat
set out = GT95_d.out

echo $data

#execute code
set conjugate = "-C10"

$EXEC -D$data $conjugate -d1 -M10000/100/75/10 -V -R1 -B-5/5 \
      -Pa-0/0/1/6/0/600/0/40 -Pb-0/0/1/6/0/600/0/40  > ./$out &

sleep 2
