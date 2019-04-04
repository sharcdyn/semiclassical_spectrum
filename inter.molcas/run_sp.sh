if [ ! -d SAVES ]; then
 mkdir SAVES
fi

if [ ! -e SAVES/sp_$1.tar.bz2 ]; then
 run_openmolcas sp
 python3 take_sp.py sp.rassi.h5
 tar cjf SAVES/sp_$1.tar.bz2 H0.dat dipole*.dat sp.out
else
 tar xjf SAVES/sp_$1.tar.bz2
fi


