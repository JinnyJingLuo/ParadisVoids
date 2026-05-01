# Declaration of variables
CC = mpicxx
EZMATH_DIR = ../ezmath/
DISLOCATION_DIR = ../dislocation_base/
FLAGS = -I$(EZMATH_DIR) -I$(DISLOCATION_DIR)
OPTIMIZATION = -O3
EXEC = sadisgen

EZMATH_SOURCES =			Point.cpp     \
					Vector.cpp	\
					GenericNode.cpp	\
					Randomizer.cpp	\
					Tools.cpp	\
					Geometry.cpp	\
					GeometricComponent.cpp	\
					Polyhedron.cpp	\
					PolyhedronEdge.cpp	\
					PolyhedronFace.cpp	\
					AxisAlignedBoundingBox.cpp	\
					CartesianOrthogonalCoordinateSystem.cpp	\
					Matrix.cpp	\
					Patch.cpp	\
					TriPatch.cpp	\
					Line.cpp	\
					Segment.cpp	\
					Plane.cpp

DISLOCATION_SOURCES = PrecipitateStructure.cpp

SOURCES = 	Chain.cpp \
			PrecipitateDislocationGenerator.cpp	\
			SADisGen.cpp


ABSOLUTE_EZMATH_SOURCES = $(foreach source,$(EZMATH_SOURCES),$(EZMATH_DIR)$(source))
EZMATH_OBJECTS = $(ABSOLUTE_EZMATH_SOURCES:.cpp=.o)
ABSOLUTE_DISLOCATION_SOURCES = $(foreach source,$(DISLOCATION_SOURCES),$(DISLOCATION_DIR)$(source))
DISLOCATION_OBJECTS = $(ABSOLUTE_DISLOCATION_SOURCES:.cpp=.o)
OBJECTS = $(SOURCES:.cpp=.o)

# Main target
$(EXEC): $(EZMATH_OBJECTS) $(DISLOCATION_OBJECTS) $(OBJECTS)
	$(CC) $(EZMATH_OBJECTS) $(DISLOCATION_OBJECTS) $(OBJECTS) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(FLAGS) $(OPTIMIZATION) $< -o $@


# To remove generated files
clean:
	rm -f $(EXEC) ./*.o
        
        


