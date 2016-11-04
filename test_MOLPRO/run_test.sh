#!/bin/bash

MOLP=molpros_2012.1
TMPDIR=$PWD

if [ ! -e ref.wf ]; then
 $MOLP -W$TMPDIR -I$TMPDIR -d$TMPDIR cas0
fi

../bin/sharc_spec.exe < spec.in > spec.out
