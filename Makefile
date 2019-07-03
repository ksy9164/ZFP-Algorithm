CC = gcc
TARGET = zfp
SOURCES = main.c zfp.c 
LDFLAGS = -std=c11
all:
	$(CC) -o  $(TARGET) $(SOURCES) $(LDFLAGS)
clean:
	rm -rf *.o $(TARGET)

