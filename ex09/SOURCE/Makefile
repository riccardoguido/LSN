CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o genetic.o
	$(CC) genetic.o random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
genetic.o : genetic.cpp genetic.h 
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe
remove:
	rm ../OUTPUT/*.*
