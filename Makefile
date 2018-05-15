include Make.inc
CFLAGS-add += -DMATHLIB_STANDALONE

SRC = $(wildcard *.c)
TAR = $(SRC:.c=.so)

all: $(TAR)

%.so: %.c
	$(CC) $(CFLAGS-add) -march=native -o $@ $^