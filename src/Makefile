FLAGS           =
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =
CLEANFILES       =

### This can be defined from outside:
# make COMP=gfortran
UNAME = ${shell uname}

COMP=gfortran
STATIC=no

F77C= $(COMP)
F90C= $(COMP)
F90OPTS=-O2
F77OPTS=-O2
LINKOPTS=-I. -llapack

OBJS=getdata.o gauss.o diag-lapack.o

MKDIR_P=mkdir -p
BIN=../bin
.PHONY: ${BIN}

all: ${BIN} spec

${BIN}:
	${MKDIR_P} ${BIN}

spec: $(OBJS) sharc_spec.o
	$(F90C) $(F90OPTS) -o ${BIN}/sharc_spec.exe sharc_spec.o $(OBJS) \
	$(LINKOPTS)

clean:
	rm -f *.o *.mod

.SUFFIXES : .f90 .f

.f90.o:
	$(F90C) $(F90OPTS) -c $<

.f.o:
	$(F77C) $(F77OPTS) -c $<

