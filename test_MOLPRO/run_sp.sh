#!/bin/bash

MOLP=molpros_2012.1
TMPDIR=$PWD

traj=` echo $1 | awk '{printf "%10.10i",$1}' `

if [ ! -e $traj.tar.bz2  ]; then
 $MOLP -W$TMPDIR -I$TMPDIR -d$TMPDIR cas
 tar cjf $traj.tar.bz2 cas.out cas.molden
 rm cas.*
fi

 tar xjf $traj.tar.bz2 

### Extracting Spin-Orbit matrix
 awk '/Lowest unperturbed energy/{E0=$NF}
     / Spin-Orbit Matrix /{pr=1}
     /Nr  State  S   SZ/{nst=0;typdat=2;npot=$NF;for (i=5;i<=NF;i++) {coord1[i-4]=$i}}
     {if (pr==1){
      if (typdat==1) {typdat=2;
       for (i=1;i<=NF;i++) {SO[typdat,nst,coord1[i]]=$i/2.194746E5}}
      if ($1==nst+1 && typdat==2) {typdat=1;nst++;
       for (i=5;i<=NF;i++) {SO[typdat,nst,coord1[i-4]]=$i/2.194746E5}}
     }}END{
      for (typdat=1;typdat<=2;typdat++){
       if (typdat==1) {fich="H0.dat"}
       if (typdat==2) {fich="H0.cmp"}
       for (i=1;i<=npot;i++){
        if (typdat==1) {SO[1,i,i]+=E0}
        for (j=1;j<=npot;j++){
         printf " %30.20e",SO[typdat,i,j] > fich
        }
        printf "\n" > fich
       }
     }}' cas.out

### Extracting SO-dipole moments (can be done also with the Spin-Free ones using the files dipole1.dat, dipole2.dat dipole3.dat)
### This can have some problems if there is degenartion, since the dipoles are asigned by energy
 awk '
     /Property matrix of the DMX operator/{dip=1}
     /Property matrix of the DMY operator/{dip=2}
     /Property matrix of the DMZ operator/{dip=3}
     /Nr Sym/{pr=1;nst=0;for (i=5;i<=NF;i++){coord1[i-4]=$i};npot=$NF}
     /Expectation/{pr=0}
     {if (pr==1) {
      if ($1==nst+1) {typdat=1;nst++;
       for (i=3;i<=NF;i++){dipo[1,dip,nst,coord1[i-2]]=$i}}
      if (typdat==1) {typdat=0;
       for (i=3;i<=NF;i++){dipo[2,dip,nst,coord1[i-2]]=$i}}
     }}END{
      for (typdat=1;typdat<=2;typdat++){
       for (dip=1;dip<=3;dip++){
        if (typdat==1) {ext="r"}
        if (typdat==2) {ext="i"}
        if (dip==1) {header="1"}
        if (dip==2) {header="2"}
        if (dip==3) {header="3"}
        fich="dip"header""ext".dat"
        for (i=1;i<=npot;i++){
         for (j=1;j<=npot;j++){
          printf " "dipo[typdat,dip,i,j] > fich
         }
         printf "\n" > fich
        }
     }}}' cas.out

 rm cas.*
