CC = gcc
CFLAGS = -Wall -Iinclude
SRC = src/main.c src/q1.c src/q2.c
OUT = programa

all: $(OUT)

$(OUT): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(OUT)

clean:
	rm -f $(OUT)
