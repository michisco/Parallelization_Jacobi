.PHONY: clean

CXX	 = g++ -std=c++20 -O3
CXXFLAGS = -pthread
SRC 	 = ./src
ALL	 = jacobi_seq jacobi_par jacobi_pinned jacobi_ff 

all: $(ALL)

jacobi_pinned: jacobi_main_pinned.cpp $(SRC)/parallelJacobi.h $(SRC)/sequentialJacobi.h $(SRC)/utimer.h
	$(CXX) $(CXXFLAGS) -I $(SRC) $< -o $@
	
jacobi_par: jacobi_main_barrier.cpp $(SRC)/parallelJacobi.h $(SRC)/sequentialJacobi.h $(SRC)/utimer.h
	$(CXX) $(CXXFLAGS) -I $(SRC) $< -o $@

jacobi_ff: jacobi_main_fastflow.cpp $(SRC)/fflowJacobi.h $(SRC)/sequentialJacobi.h $(SRC)/utimer.h
	$(CXX) $(CXXFLAGS) -I $(SRC) $< -o $@

jacobi_seq: jacobi_main_sequential.cpp $(SRC)/sequentialJacobi.h $(SRC)/utimer.h
	$(CXX) $(CXXFLAGS) -I $(SRC) $< -o $@

clean:
	-rm $(ALL)
	-rm *.o
