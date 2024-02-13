LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

do: 
	make Directionality_ToyMC
	make Ordered_hits

Directionality_ToyMC: Directionality_ToyMC.o Functions.o
	g++ Directionality_ToyMC.o Functions.o -o Directionality_ToyMC -lMinuit ${INCS} ${LIBS}
Ordered_hits: Ordered_hits.o Functions.o
	g++ Ordered_hits.o Functions.o -o Ordered_hits -lMinuit ${INCS} ${LIBS}

Ordered_hits.o: Ordered_hits.cxx 
	g++ -c Ordered_hits.cxx -o Ordered_hits.o -lMinuit ${INCS} ${LIBS}
Directionality_ToyMC.o: Directionality_ToyMC.cxx 
	g++ -c Directionality_ToyMC.cxx -o Directionality_ToyMC.o -lMinuit ${INCS} ${LIBS}
Functions.o: Functions.cxx Functions.h 
	g++ -c Functions.cxx -o Functions.o ${INCS}

clean:
	rm *.o Output.*
