CXXFLAGS = -std=c++11 -Wall -g

all : projection optimization multiOpt

projection : projection.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) -O2 -larmadillo -llapack -lblas

optimization : optimization_test.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) -O2 -larmadillo -llapack -lblas

multiOpt : multi_test.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) -O2 -larmadillo -llapack -lblas

cellTest : cell_testing.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) -O2 -larmadillo -llapack -lblas

testProj : projection
	./projection
	python displayTest.py

testOpt : optimization
	./optimization
	python opt_display.py

testMulti : multiOpt
	./multiOpt
	python multiOpt_display.py

clean :
	rm -f *.o *~

.PHONY : all clean
