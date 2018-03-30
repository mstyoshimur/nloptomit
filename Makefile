#Makefile
nloptomit: nloptomit.o
	g++  nloptomit.o -L${CCP4}/lib -lclipper-core   -lclipper-ccp4 -lclipper-cif -lclipper-minimol -lclipper-cns -lclipper-contrib   -lclipper-mmdb -lclipper-phs -lfftw -lmmdb2 -lccp4c  -lgsl -lgslcblas -lnlopt  -lm -o nloptomit 

nloptomit.o: nloptomit.cpp
	g++ -c nloptomit.cpp -I${CCP4src}/ccp4-7.0-src/checkout/clipper -I${CCP4}/include 