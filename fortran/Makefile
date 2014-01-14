#2DMHD Makefile
#Define Objects

OBJ = obj
SRC = src
MOD = mod
PLOTS = plots
DATA = data
EXEC = xmhd2d

OBJECTS = ${OBJ}/global.o\
					${OBJ}/vector_operations.o\
					${OBJ}/interpolation.o\
					${OBJ}/gen_grid.o\
					${OBJ}/get_Va.o\
					${OBJ}/get_facts.o\
					${OBJ}/get_basisfn.o\
					${OBJ}/main.o


#Define Compiler

#############ifort############
COMP = ifort

#Define Switches
#SWITCH = -O2 -module ${MOD}
SWITCH = -g -debug-parameters -module ${MOD} -check all

#Define Libraries
LIB_DIR = -L lib/ifort
LIB = -llapack -lblas

##############################


#############gfortran############
#COMP = gfortran

#Define Switches
#SWITCH = -O2 -J ${MOD}
#SWITCH = -g -J ${MOD}

#Define Libraries

#LIB_DIR = -L lib/gfortran
#LIB = -llapack -lblas

##############################


#Makefile
all: DIR $(OBJECTS)
	$(COMP) $(SWITCH) ${OBJECTS} -o ${EXEC} $(LIB_DIR) $(LIB)

${OBJ}/global.o: ${SRC}/global.f90
	$(COMP) -c $(SWITCH) ${SRC}/global.f90 -o ${OBJ}/global.o

${OBJ}/main.o: ${OBJ}/global.o ${SRC}/main.f90 ${OBJ}/gen_grid.o ${OBJ}/get_Va.o ${OBJ}/get_facts.o
	$(COMP) -c $(SWITCH) ${SRC}/main.f90 -o ${OBJ}/main.o

${OBJ}/gen_grid.o: ${OBJ}/global.o ${SRC}/gen_grid.f90
	$(COMP) -c $(SWITCH) ${SRC}/gen_grid.f90 -o ${OBJ}/gen_grid.o

${OBJ}/get_Va.o: ${OBJ}/global.o ${SRC}/get_Va.f90 ${OBJ}/interpolation.o
	$(COMP) -c $(SWITCH) ${SRC}/get_Va.f90 -o ${OBJ}/get_Va.o

${OBJ}/interpolation.o: ${SRC}/interpolation.f90
	$(COMP) -c $(SWITCH) ${SRC}/interpolation.f90 -o ${OBJ}/interpolation.o

${OBJ}/get_facts.o: ${SRC}/get_facts.f90 ${OBJ}/global.o
		$(COMP) -c $(SWITCH) ${SRC}/get_facts.f90 -o ${OBJ}/get_facts.o

${OBJ}/get_basisfn.o: ${SRC}/get_basisfn.f90 ${OBJ}/global.o
		$(COMP) -c $(SWITCH) ${SRC}/get_basisfn.f90 -o ${OBJ}/get_basisfn.o

${OBJ}/vector_operations.o: ${SRC}/vector_operations.f90
	$(COMP) -c $(SWITCH) ${SRC}/vector_operations.f90 -o ${OBJ}/vector_operations.o

DIR:
	mkdir -p $(OBJ) $(DATA) $(PLOTS) $(MOD)
clean:
	rm -rf ${OBJ}

#Cleaning Everything
cleanest: clean
	rm -rf ${EXEC}
	rm -rf ${DATA}
	rm -rf ${PLOTS}
	rm -rf ${MOD}

# End of the makefile

