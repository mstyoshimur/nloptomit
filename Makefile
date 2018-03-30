#Makefile
cnlpomit: cnlpomit.o
	g++  cnlpomit.o -L/work/ccp4/ccp4-7.0-src/ccp4-7.0/lib -lclipper-core   -lclipper-ccp4 -lclipper-cif -lclipper-minimol -lclipper-cns -lclipper-contrib   -lclipper-mmdb -lclipper-phs -lfftw -lmmdb2 -lccp4c  -lgsl -lgslcblas -lnlopt  -lm -o cnlpomit 

cnlpomit.o: cnlpomit.cpp
	g++ -c cnlpomit.cpp -I/work/ccp4/ccp4-7.0-src/checkout/clipper -I/work/ccp4/ccp4-7.0-src/ccp4-7.0/include 