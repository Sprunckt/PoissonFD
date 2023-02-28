CFLAGS = -O3 -Wall -Wextra
CC = g++
all: main

main: main.o linalg.o model.o mesh.o
	$(CC) $(CFLAGS)  main.o linalg.o model.o mesh.o -o main

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

model.o: model.cpp
	$(CC) $(CFLAGS)  -c model.cpp
	
mesh.o: mesh.cpp
	$(CC) $(CFLAGS)  -c mesh.cpp

linalg.o: linalg.cpp
	$(CC) $(CFLAGS)  -c linalg.cpp


clean:
	rm -f *.o main linalg model 