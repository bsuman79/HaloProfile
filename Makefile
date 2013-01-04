CC = g++
OFILES=mpfit.o curvefit.o concentration.o solve_c.o fsigma.o sphericalcollapse.o bin_conc.o driver_conc.o
LIBFILE=libmpfit.a

FLAGS_INCL = -I/usr/local/include
LIB_PATH = -L/usr/local/lib -lgsl -lgslcblas -lmpfit -lm


image  = driver_conc
#curvefit

all: $(LIBFILE) $(image)


clean:
	rm -f $(OFILES) $(image) $(LIBFILE) driver_mf driver_conc

mpfit.o: mpfit.c mpfit.h
	$(CC) -c -o $@ $< $(CFLAGS)

sphericalsollapse.o: sphericalcollapse.c
	$(CC) -c -o $@ $< $(CFLAGS)

fsigma.o: fsigma.cxx
	$(CC) -c -o $@ $< $(CFLAGS)

curvefit.o: curvefit.c 
	$(CC) -c -o $@ $< $(CFLAGS)

solve_c.o: solve_c.c
	$(CC) -c -o $@ $< $(CFLAGS)

concentration.o: concentration.cxx
	$(CC) -c -o $@ $< $(CFLAGS)

bin_conc.o: bin_conc.cxx
	$(CC) -c -o $@ $< $(CFLAGS)

driver_conc.o: driver_conc.cxx
	$(CC) -c -o $@ $< $(CFLAGS)

$(LIBFILE): $(OFILES)
	$(AR) r $@ $(OFILES)

$(image): $(image).cxx libmpfit.a
	g++  -g $(FLAGS_INCL) $(image).o solve_c.o curvefit.o fsigma.o sphericalcollapse.o mpfit.o bin_conc.o concentration.o $(CFLAGS) -o $(image) -L.  $(LIB_PATH)

