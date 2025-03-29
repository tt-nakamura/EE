NTL = -lntl -lgmp -L/usr/local/lib
OBJ = EE.o EEFactoring.o ZZlib.o ZZFactoring.o mpqs.o rho.o

example: example.o CubRootMod.o $(OBJ)
	g++ example.o CubRootMod.o $(OBJ) $(NTL)
fig1: fig1.o $(OBJ)
	g++ fig1.o $(OBJ) $(NTL)
