#2DMHD Makefile
#Define Objects

OBJ = obj
SRC = src
MOD = mod
EXEC = xmhd2d

OBJECTS = ${OBJ}/global.o\
			${OBJ}/main.o\
			${OBJ}/gen_grid.o			
#Define Compiler

#############ifort############
COMP = ifort
	
#Define Switches
#SWITCH = -O2 -module ${MOD}
SWITCH = -g -debug-parameters -module ${MOD}

#Define Libraries
LIB_DIR = -L lib/ifort
LIB = -llapack -lblas

##############################


#############gfortran############
#COMP = gfortran
	
#Define Switches
#SWITCH = -O2 -J${MOD}
#SWITCH = -g -J${MOD}

#Define Libraries

#LIB_DIR = -L lib/gfortran
#LIB = -llapack -lblas

##############################


#Makefile
all: $(OBJECTS)
	$(COMP) $(SWITCH) ${OBJECTS} -o ${EXEC} $(LIB_DIR) $(LIB)
	
${OBJ}/global.o: ${SRC}/global.f90
	$(COMP) -c $(SWITCH) ${SRC}/global.f90 -o ${OBJ}/global.o
	
${OBJ}/main.o: ${MOD}/global.mod ${SRC}/main.f90
	$(COMP) -c $(SWITCH) ${SRC}/main.f90 -o ${OBJ}/main.o
${OBJ}/gen_grid.o: ${MOD}/global.mod ${SRC}/gen_grid.f90
	$(COMP) -c $(SWITCH) ${SRC}/gen_grid.f90 -o ${OBJ}/gen_grid.o

clean:
	rm ${OBJ}/*.o
	rm ${MOD}/*.mod
#Cleaning Everything
cleanest: clean
	rm ${EXEC}
# End of the makefile

