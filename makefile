CXXFLAGS =  -O3  -Wall -march=native -lm -I.

SU2_adjoint: ranlxd.o geometry.o update.o APE_smearing.o observables.o archive.o main_measure.o measure.o main_SU2_adjoint.o
	g++ $(CXXFLAGS) ranlxd.o geometry.o update.o APE_smearing.o observables.o archive.o main_SU2_adjoint.o -o SU2_adjoint -lm
	g++ $(CXXFLAGS) ranlxd.o geometry.o update.o APE_smearing.o observables.o archive.o measure.o -o measure -lm
	g++ $(CXXFLAGS) ranlxd.o geometry.o update.o APE_smearing.o observables.o archive.o main_measure.o -o main_measure -lm
ranlxd.o: ranlxd.c
	g++ $(CXXFLAGS) -c ranlxd.c
geometry.o: geometry.cpp
	g++ $(CXXFLAGS) -c geometry.cpp
archive.o: archive.cpp
	g++ $(CXXFLAGS) -c archive.cpp
update.o: update.cpp
	g++ $(CXXFLAGS) -c update.cpp
observables.o: update.o
	g++ $(CXXFLAGS) -c observables.cpp
APE_smearing.o: APE_smearing.cpp
	g++ $(CXXFLAGS) -c APE_smearing.cpp
main_SU2_adjoint.o: main_SU2_adjoint.cpp
	g++ $(CXXFLAGS) -c main_SU2_adjoint.cpp
measure.o: measure.cpp
	g++ $(CXXFLAGS) -c measure.cpp
main_measure.o: main_measure.cpp
	g++ $(CXXFLAGS) -c main_measure.cpp
prepare:
	mkdir conf
	mkdir adj_conf
	mkdir O1minus_output_files
	mkdir O0plus_output_files
	mkdir log_files
	mkdir obs
clear:
	rm -rf *.o
