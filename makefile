OBJS = simutrans.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

simutrans : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o simutrans

clean:
	\rm *.o simutrans
