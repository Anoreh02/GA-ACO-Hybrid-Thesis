CC=gcc
CFLAGS=-O2 -std=c11 -Iinclude
SRC=$(wildcard src/*.c)
BIN=bin/gaacosp

all: $(BIN)

$(BIN): $(SRC)
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $(SRC)

clean:
	rm -rf bin
