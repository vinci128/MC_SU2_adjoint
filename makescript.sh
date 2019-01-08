mkdir conf
mkdir adj_conf
g++ -Wall -g -c ranlxd.c
g++ -Wall -g -c geometry.c
g++ -Wall -g -c archive.c
g++ -Wall -g -c update.c
g++ -Wall -g -c observables.c
g++ -Wall -g -c APE_smearing.c
g++ -Wall -g -c main_SU2_adjoint.c
g++ -Wall -g ranlxd.o geometry.o update.o APE_smearing.o observables.o archive.o main_SU2_adjoint.o -o SU2_adjoint -lm
