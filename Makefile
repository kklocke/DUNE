CXXFLAGS = -std=c++11 -Wall

all : projection optimization

projection : projection.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) -O2 -larmadillo -llapack -lblas

optimization : optimization_test.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) -O2 -larmadillo -llapack -lblas

clean :
	rm -f *.o *~

.PHONY : all clean
