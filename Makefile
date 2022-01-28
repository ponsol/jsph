CPP=g++
FC=gfortran
FCFLAGS= -O1

OBJS= eos.o part.o kernel.o problem.o  artvis.o  geom.o sphsub.o jsph.o

all: jsph

jsph: ${OBJS}
	$(FC) $(FCFLAGS) -o jsph ${OBJS}

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
	

%.mod: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -rf ${OBJS} jsph *.mod pall pcut
