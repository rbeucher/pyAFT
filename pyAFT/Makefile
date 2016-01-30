
CC = gcc
CFLAGS = -c -fpic
LDFLAGS = -shared

SRC_DIR = ./src
SRC_FILES = $(SRC_DIR)/Ketcham.c
INCLUDE = $(SRC_DIR)/ketch.h

DST_DIR = ./lib


OBJECTS = $(patsubst %.c,%.o,$(SRC_FILES)) 

.SUFFIXES: .o .c .h

.c.o:
	$(CC) $(CFLAGS) -I$(SRC_DIR) $*.c
	mv *.o $(SRC_DIR)/.

ketcham: $(OBJECTS) 
	$(CC) $(LDFLAGS) -o $(DST_DIR)/libketcham.so $(OBJECTS)

clean:
	rm $(SRC_DIR)/*.o
	rm $(DST_DIR)/*.so
