CXXFLAGS = -std=c++11 -Wall

all : projection

projection : projection.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean :
	rm -f *.o *~

.PHONY : all clean
