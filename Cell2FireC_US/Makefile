EIGENDIR     = /home/minho/research/fire/c2f_codes/Eigen34
CC = g++
MPCC = g++
OPENMP = -openmp
CFLAGS = -std=c++11 -O3 -I$(EIGENDIR)
LIBS = -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -lm -lpthread -ldl 

TARGETS = Cell2Fire

all:	$(TARGETS)

Cell2Fire: Cell2Fire.o CellsFBP.o fuelmodelBP.o SpottingFBP.o ReadCSV.o ReadArgs.o WriteCSV.o Ellipse.o
	$(CC) -o $@ $(LIBS) -fopenmp Cell2Fire.o CellsFBP.o fuelmodelBP.o SpottingFBP.o ReadCSV.o ReadArgs.o WriteCSV.o Ellipse.o

Cell2Fire.o: Cell2Fire.cpp CellsFBP.o fuelmodelBP.o SpottingFBP.o ReadCSV.o ReadArgs.o WriteCSV.o
	$(CC) -c $(CFLAGS) -fopenmp Cell2Fire.cpp

SpottingFBP.o: SpottingFBP.cpp SpottingFBP.h CellsFBP.h
	$(CC) -c $(CFLAGS) SpottingFBP.cpp CellsFBP.h

CellsFBP.o: CellsFBP.cpp CellsFBP.h fuelmodelBP.o
	$(CC) -c $(CFLAGS) CellsFBP.cpp

fuelmodelBP.o: fuelmodelBP.cpp fuelmodelBP.h
	$(CC) -c $(CFLAGS) fuelmodelBP.cpp

ReadCSV.o: ReadCSV.cpp ReadCSV.h fuelmodelBP.o
	$(CC) -c $(CFLAGS) ReadCSV.cpp

ReadArgs.o: ReadArgs.cpp ReadArgs.h
	$(CC) -c $(CFLAGS) ReadArgs.cpp

# Lightning.o: Lightning.cpp Lightning.h
# 	$(CC) -c $(CFLAGS) Lightning.cpp

Forest.o: Forest.cpp Forest.h
	$(CC) -c $(CFLAGS) Forest.cpp

WriteCSV.o: WriteCSV.cpp WriteCSV.h
	$(CC) -c $(CFLAGS) WriteCSV.cpp

Ellipse.o: Ellipse.cpp Ellipse.h 
	$(CC) -c $(LIBS) $(CFLAGS) Ellipse.cpp 


clean:
	rm  ReadArgs.o ReadCSV.o fuelmodelBP.o Cell2Fire.o CellsFBP.o Cell2Fire SpottingFBP.o Forest.o WriteCSV.o Ellipse.o *.gch