# Example makefile, please adjust for your system.

CC = clang
OPT = -O2 -std=c99 -g -Wall
LINK = -L/usr/local/lib -L/usr/lib -L/opt/local/lib
INCLUDE = -I. -I/usr/include -I/opt/local/include
LIB = -lm -lgsl -lgslcblas
FLAGS = ${OPT} ${DEBUG} ${INCLUDE}

OBJECTS = dd.o input.o psvindex.o jobfile.o output.o forces.o rotations.o initialvalues.o geomelem.o impact.o odeint.o rkdopr853qs.o rkdopr853.o


bin: ${OBJECTS}
	${CC} ${OBJECTS} ${FLAGS} ${LINK} ${LIB} -o dd.x

$(OBJECTS): %.o: %.c
	${CC} ${FLAGS} -c -o $@ $<

clean:
	rm *.o
