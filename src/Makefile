CC = gcc
CFLAGS = -Wall -O3
LFLAGS = -lexpat -lm

EXEC = cpg
MAIN = cpg.c
OBJS = 

all: $(EXEC)

.cc.o:
	$(CC) $(CFLAGS) -c $<
	
$(EXEC): $(MAIN) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

release: $(EXEC) clean
	strip $(EXEC)

clean:
	rm -f $(OBJS) core

distclean: clean
	rm -f $(EXEC)
