INCLUDE =
CC = mpicc
CFLAGS = -O2 -lm

%.o:	%.c
	$(CC) -c $^ $(INCLUDE) $(CFLAGS)

mvm:	mvm.o mvm-driver.o
	$(CC) -o $@ $^ $(CFLAGS)

tdmvm:	tdmvm.o tdmvm-driver.o
	$(CC) -o $@ $^ $(CFLAGS)
