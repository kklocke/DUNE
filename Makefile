CXXFLAGS = -std=c++14 -Wall

all : projection

projection : projection.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean :
	rm -f test-maze *.o *~

.PHONY : all clean
