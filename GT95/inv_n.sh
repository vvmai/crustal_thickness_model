#!/bin/csh -f
#
# inversion of Gleason & Tullis (1994) data on dry dislocation creep

#set up execution
set EXEC = labmc

set data = GT95_V_corr.dat
set out = GT95_V_corr.out

echo $data

#execute code
#set conjugate = "-C10"
set conjugate = ""

$EXEC -D$data $conjugate -M1000000/100/75/10 -V -R1 -B-5/5 \
      -Pb-0/0/1/6/0/600/0/40  > ./$out &

sleep 2

end