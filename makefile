CXX=g++ -std=c++17
OPT=-O3

CSP: a_FHL_rec.o -lboost_system -lboost_thread
	$(CXX) -g a_FHL_rec.o -lboost_system -lboost_thread

a_FHL_rec.o:a_FHL_rec.cpp
	$(CXX) -g -c $(OPT) a_FHL_rec.cpp 

clean:
	rm *.o
	rm a.out
