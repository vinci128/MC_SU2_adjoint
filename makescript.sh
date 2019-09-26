mkdir conf
mkdir adj_conf
mkdir O1minus_output_files
mkdir O0plus_output_files
mkdir log_files
mkdir obs
g++ -Wall -g -c ranlxd.c
g++ -Wall -g -c geometry.cpp
g++ -Wall -g -c archive.cpp
g++ -Wall -g -c update.cpp
g++ -Wall -g -c observables.cpp
g++ -Wall -g -c APE_smearing.cpp
g++ -Wall -g -c main_SU2_adjoint.cpp
g++ -Wall -g -c measure.cpp
g++ -Wall -g -c main_measure.cpp
g++ -Wall -g ranlxd.o geometry.o update.o APE_smearing.o observables.o archive.o main_SU2_adjoint.o -o SU2_adjoint -lm
g++ -Wall -g ranlxd.o geometry.o update.o APE_smearing.o observables.o archive.o measure.o -o measure -lm
g++ -Wall -g ranlxd.o geometry.o update.o APE_smearing.o observables.o archive.o main_measure.o -o main_measure -lm
