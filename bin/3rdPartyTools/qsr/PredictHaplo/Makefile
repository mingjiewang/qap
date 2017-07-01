
CFLAGS  =	-O3 -DSCYTHE_LAPACK -DSCYTHE_DEBUG=0 -DHAVE_TRUNC -fPIC 

SCYTHE  =	NEWSCYTHE/include/

EXOBJS	=	PredictHaplo_externAlign.o


all:	PredictHaplo-Paired	

clean:
	rm $(EXOBJS) PredictHaplo_externAlign

PredictHaplo-Paired:	$(EXOBJS)
	g++ $(CFLAGS)  -o $@  $(EXOBJS)  -lblas -llapack
PredictHaplo_externAlign.o:	PredictHaplo_externAlign.cpp 
	g++   $(CFLAGS)  -I$(SCYTHE) -c -o $@  $<
