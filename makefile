LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

do: Directionality_ToyMC.cxx Ordered_hits.cxx
	make Directionality_ToyMC
	make Ordered_hits

Directionality_ToyMC: Directionality_ToyMC.cxx
	g++ Directionality_ToyMC.cxx -o Directionality_ToyMC -lMinuit ${INCS} ${LIBS}
Ordered_hits: Ordered_hits.cxx
	g++ Ordered_hits.cxx -o Ordered_hits -lMinuit ${INCS} ${LIBS}

clean:
	rm *.o Output.*
