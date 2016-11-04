subroutine getdata(npot,Epot,dipole,ntraj,nat,geom,sat,U)
implicit none
integer npot, ntraj,nat
complex*16 dipole(3,npot,npot)
real*8 Epot(npot),geom(nat,3)
character*2 sat(nat)

complex*16 U(npot,npot)

real*8 Hr(npot,npot),Hi(npot,npot),Ur(npot,npot),Ui(npot,npot)
real*8 Ar(npot,npot),Ai(npot,npot)
logical ex0,ex
integer i,j,k
real*8 fv1(npot),fv2(npot),fm1(2,npot)
integer ierr

logical ex1(3),ex2(3)
character*100 exe
character*30 fich1(3),fich2(3)

open(1,file="geom.xyz")
write(1,*) nat
write(1,*)
do i=1,nat
 write(1,"(A2,3(x,F20.10))") sat(i),(geom(i,j)*.5292,j=1,3)
enddo
close(1)

write(exe,"(A,I6)") "bash run_sp.sh ",ntraj
call system(exe)

inquire(file="H0.dat",exist=ex0)
if (ex0) then
 open(1,file="H0.dat")
 do i=1,npot
  read(1,*) (Hr(i,j),j=1,npot)
 enddo
 close(1)

 inquire(file="H0.cmp",exist=ex)
 if (ex) then
  open(1,file="H0.cmp")
  do i=1,npot
   read(1,*) (Hi(i,j),j=1,npot)
  enddo
  close(1)
 else
  do i=1,npot
   do j=1,npot
    Hi(i,j)=0.d0
   enddo
  enddo
  if (ntraj.eq.0) write(6,*) "H0.cmp does not exist for equilibrium geometry"
 endif

! Calculating the U matrix for H0
 Ar=Hr
 Ai=Hi
 U=dcmplx(Ar,Ai)
 call cdiag(U,npot,Epot)
 call corder(U,npot,Epot,.false.)
 Ur=dreal(U)
 Ui=dimag(U)
 if (ntraj.eq.0) then
  write(6,*) "CAS energies"
  write(6,*) Epot(:)
 endif

 do k=1,3
  do i=1,npot
   do j=1,npot
    dipole(k,i,j)=dcmplx(0.d0,0.d0)
   enddo
  enddo
 enddo

 fich1(1)="dip1r.dat"
 fich1(2)="dip2r.dat"
 fich1(3)="dip3r.dat"
 fich2(1)="dip1i.dat"
 fich2(2)="dip2i.dat"
 fich2(3)="dip3i.dat"
 ex=.false.
 do k=1,3
  inquire(file=fich1(k),exist=ex1(k))
  inquire(file=fich2(k),exist=ex2(k))
  if (ex1(k).or.ex2(k)) ex=.true.
 enddo
 if (ex) then
  if (ntraj.eq.0) then
   write(6,*) "Trying adiabatic dipoles"
   do i=1,3
    write(6,*) fich1(i),ex1(i),fich2(i),ex2(i)
   enddo
  endif
  do k=1,3
   if (ex1(k)) open(1,file=fich1(k))
   if (ex2(k)) open(2,file=fich2(k))
   do i=1,npot
    do j=1,npot
     Ar(i,j)=0.d0
     Ai(i,j)=0.d0
    enddo
    if (ex1(k)) read(1,*) (Ar(i,j),j=1,npot)
    if (ex2(k)) read(2,*) (Ai(i,j),j=1,npot)
   enddo
   call diabatize (npot,Ar,Ai,Ur,Ui)
   do i=1,npot
    do j=1,npot
     dipole(k,i,j)=dcmplx(Ar(i,j),Ai(i,j))
    enddo
   enddo
   if (ex1(k)) close(1)
   if (ex2(k)) close(2)
  enddo

 else

   ! Trying with diabatic dipoles (SHARC propagation style)
  fich1(1)="dipole1.dat"
  fich1(2)="dipole2.dat"
  fich1(3)="dipole3.dat"
  fich2(1)="dipole1.cmp"
  fich2(2)="dipole2.cmp"
  fich2(3)="dipole3.cmp"
  do k=1,3
   inquire(file=fich1(k),exist=ex1(k))
   inquire(file=fich2(k),exist=ex2(k))
  enddo
  if (ntraj.eq.0) then
   write(6,*) "Trying diabatic dipoles"
   do i=1,3
    write(6,*) fich1(i),ex1(i),fich2(i),ex2(i)
   enddo
  endif
  do k=1,3
   if (ex1(k)) open(1,file=fich1(k))
   if (ex2(k)) open(2,file=fich2(k))
   do i=1,npot
    do j=1,npot
     Ar(i,j)=0.d0
     Ai(i,j)=0.d0
    enddo
    if (ex1(k)) read(1,*) (Ar(i,j),j=1,npot)
    if (ex2(k)) read(2,*) (Ai(i,j),j=1,npot)
    do j=1,npot
     dipole(k,i,j)=dcmplx(Ar(i,j),Ai(i,j))
    enddo
   enddo
   if (ex1(k)) close(1)
   if (ex2(k)) close(2)
  enddo
 endif

 ! Looking for real Hamiltonian
 inquire(file="H1.dat",exist=ex1(1))
 inquire(file="H1.cmp",exist=ex2(1))
 if (ex1(1).or.ex2(1)) then
 if (ex1(1)) open(1,file="H1.dat")
  if (ex2(1)) open(2,file="H1.cmp")
  do i=1,npot
   do j=1,npot
    Ar(i,j)=0.d0
    Ai(i,j)=0.d0
   enddo
   if (ex1(1)) read(1,*) (Ar(i,j),j=1,npot)
   if (ex2(1)) read(2,*) (Ai(i,j),j=1,npot)
   do j=1,npot
    if (i.eq.j) then
     Hr(i,j)=Ar(i,j)
     Hi(i,j)=Ai(i,j)
    else
     Hr(i,j)=Hr(i,j)+Ar(i,j)
     Hi(i,j)=Hi(i,j)+Ai(i,j)
    endif
   enddo
  enddo
  if(ex1(1)) close(1)
  if(ex2(1)) close(2)
 endif

 Ar=Hr
 Ai=Hi
 U=dcmplx(Ar,Ai)
 call cdiag(U,npot,Epot)
 call corder(U,npot,Epot,.false.)
 Ur=dreal(U)
 Ui=dimag(U)
 call adiabatize(npot,Hr,Hi,Ur,Ui)
 do i=1,npot
  do j=1,npot
   U(i,j)=Ur(i,j)+(0.,1.)*Ui(i,j)
  enddo
 enddo
 if (ntraj.eq.0) then
  write(6,*) "Real energies"
  write(6,*) Epot(:)
 endif

 !Checking adiabatization
! if (ntraj.eq.0) then
!  do i=1,npot
!   do j=1,npot
!    if (i.eq.j) write(90,*) "Diag ",Epot(i)
!    write(90,*) i,j,Hr(i,j),Hi(i,j)
!   enddo
!  enddo
! endif

 ! DIPOLES

 do k=1,3
  do i=1,npot
   do j=1,npot
    Ar(i,j)=dreal(dipole(k,i,j))
    Ai(i,j)=dimag(dipole(k,i,j))
   enddo
  enddo
  call adiabatize(npot,Ar,Ai,Ur,Ui)
  do i=1,npot
   do j=1,npot
    dipole(k,i,j)=dcmplx(Ar(i,j),Ai(i,j))
   enddo
  enddo
 enddo

else
 if (ntraj.eq.0) then
  write(6,*) "Must exist H0.dat for the equilibrium geometry"
  stop
 else
  write(6,*) "Skipping trajectory ",ntraj
 endif
endif


end subroutine

subroutine adiabatize(npot,Ar,Ai,Ur,Ui)
implicit none
integer npot
real*8 Ar(npot,npot),Ai(npot,npot)
real*8 Ur(npot,npot),Ui(npot,npot)
complex*16 A(npot,npot),B(npot,npot),U(npot,npot)

integer i,j,k

do i=1,npot
 do j=1,npot
  A(i,j)=dcmplx(Ar(i,j),Ai(i,j))
  U(i,j)=dcmplx(Ur(i,j),Ui(i,j))
 enddo
enddo

do i=1,npot
 do j=1,npot
  B(i,j)=dcmplx(0.d0,0.d0)
  do k=1,npot
   B(i,j)=B(i,j)+A(i,k)*U(k,j)
  enddo
 enddo
enddo
do i=1,npot
 do j=1,npot
  A(i,j)=dcmplx(0.d0,0.d0)
  do k=1,npot
   A(i,j)=A(i,j)+dconjg(U(k,i))*B(k,j)
  enddo
  Ar(i,j)=dreal(A(i,j))
  Ai(i,j)=dimag(A(i,j))
 enddo
enddo

end subroutine

subroutine diabatize(npot,Ar,Ai,Ur,Ui)
implicit none
integer npot
real*8 Ar(npot,npot),Ai(npot,npot)
real*8 Ur(npot,npot),Ui(npot,npot)
complex*16 A(npot,npot),B(npot,npot),U(npot,npot)

integer i,j,k

do i=1,npot
 do j=1,npot
  A(i,j)=dcmplx(Ar(i,j),Ai(i,j))
  U(i,j)=dcmplx(Ur(i,j),Ui(i,j))
 enddo
enddo

do i=1,npot
 do j=1,npot
  B(i,j)=dcmplx(0.d0,0.d0)
  do k=1,npot
   B(i,j)=B(i,j)+A(i,k)*dconjg(U(j,k))
  enddo
 enddo
enddo
do i=1,npot
 do j=1,npot
  A(i,j)=dcmplx(0.d0,0.d0)
  do k=1,npot
   A(i,j)=A(i,j)+U(i,k)*B(k,j)
  enddo
  Ar(i,j)=dreal(A(i,j))
  Ai(i,j)=dimag(A(i,j))
 enddo
enddo

end subroutine
