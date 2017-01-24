CC = clang
CFLAGS = -Wall -O2 -g
DEPS = print_functions.h alignment.h
OBJ = print_functions.o alignment.o main.o

%.o: %.c $(DEPS)
		$(CC) $(CFLAGS) -c -o $@ $<

main: $(OBJ)
		$(CC) $(CFLAGS) -o $@ $^

clean:
	@rm $(OBJ)
	@rm main
