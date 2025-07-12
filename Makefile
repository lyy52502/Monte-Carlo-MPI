CC = mpicc
CFLAGS = -Wall -O3 -g
BINS = malaria_sim

all: $(BINS)

malaria_sim: malaria_sim.c
	$(CC) $(CFLAGS) -o $@ malaria_sim.c

clean:
	$(RM) $(BINS)

