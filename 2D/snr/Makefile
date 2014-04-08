
.PHONY: all clean

EXES = snr

SRC = $(wildcard *.cpp)
HDR = $(wildcard *.h)
OBJ = $(SRC:.cpp=.o)

LIBS = -lm

all: $(EXES)

snr: $(OBJ)
	$(CXX) $(CXXFLAGS) -o snr $(OBJ) $(LIBS)
	
%.o: %.cpp $(HDR)
	$(CXX) $(CXXFLAGS) -c $<
	
clean:
	rm -f *.o $(EXES) ../bin/snr
