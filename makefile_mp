CC = gcc
CFLAGS = -Wall -Wextra -fopenmp
LDFLAGS = -fopenmp
EXEC = main

all: $(EXEC)

$(EXEC): main.o
	$(CC) $(LDFLAGS) -o $@ $^

main.o: main.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(EXEC) *.o
