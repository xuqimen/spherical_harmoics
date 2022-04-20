# MKLROOT = /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl

# if need to run with DEBUG mode, add -DDBUG
#CPPFLAGS = -m64 -I${MKLROOT}/include -I include/ -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -ldl -lrt -O3 -Wall -g

# if run locally, use the following (add -DDBUG to turn on debug mode)
# CPPFLAGS = -m64 -I${MKLROOT}/include -I include/ -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lrt -O3 -Wall -g

# CPPFLAGS = -m64 -I${MKLROOT}/include -I include/ -L${MKLROOT}/lib/intel64 

CFLAGS = -std=gnu99 -O3 -Wall -g

# CPPFLAGS = -m64 -I${MKLROOT}/include

## load directories
# LDFLAGS = -L${MKLROOT}/lib/intel64
## load MKL
# LDLIBS = -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lrt

LDLIBS = -lm

OBJSC = spherical_harmonics.o tools.o test.o

override CC=gcc
#override CC=mpicc

LIBBASE = test_Ylm

all: ${LIBBASE}

#mytest: $(OBJSC)
#	$(CC) $(CPPFLAGS) -o $@ $^

${LIBBASE}: $(OBJSC)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(LIBBASE) $^ $(LDLIBS)

#%.o: %.c
#	$(CC) $(CPPFLAGS) -c $<

.PHONY: clean
clean:
	rm -f  *.o ${LIBBASE}
