# Declaration of variables
CC = g++
EZMATH_DIR = ../ezmath/
FLAGS = -I$(EZMATH_DIR) -DPARALLEL
OPTIMIZATION = -O3
EXEC = fmm

EZMATH_SOURCES = 		 Point.cpp \
				 Vector.cpp \
				 Matrix.cpp \
				 Randomizer.cpp \
				 Tools.cpp \
				 Geometry.cpp \
				 CartesianOrthogonalCoordinateSystem.cpp \
				 AxisAlignedBoundingBox.cpp \

SOURCES = 	FMMTreeNode.cpp \
			FMM.cpp

ABSOLUTE_EZMATH_SOURCES = $(foreach source,$(EZMATH_SOURCES),$(EZMATH_DIR)$(source))
EZMATH_OBJECTS = $(ABSOLUTE_EZMATH_SOURCES:.cpp=.o)
OBJECTS = $(SOURCES:.cpp=.o)

# Main target
$(EXEC): $(EZMATH_OBJECTS) $(OBJECTS)
	$(CC) $(EZMATH_OBJECTS) $(OBJECTS) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(FLAGS) $(OPTIMIZATION) $< -o $@


# To remove generated files
clean:
	rm -f $(EXEC) ./*.o



