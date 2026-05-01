# Declaration of variables
CC = mpicxx
EZMATH_DIR = ../ezmath/
DISLOCATION_DIR = ../dislocation_base/
FLAGS = -I$(EZMATH_DIR) -I$(DISLOCATION_DIR)
OPTIMIZATION = -O3
EXEC = surfdef

EZMATH_SOURCES = AxisAlignedBoundingBox.cpp     \
        MathServices.cpp                        \
        Matrix.cpp                              \
        CartesianOrthogonalCoordinateSystem.cpp \
        Patch.cpp                               \
        Plane.cpp                               \
        Point.cpp                               \
        QuadPatch.cpp                           \
        Randomizer.cpp                          \
        Segment.cpp                             \
        GenericNode.cpp                         \
        Tools.cpp                               \
        TriPatch.cpp                            \
        Line.cpp                                \
        Vector.cpp                              \
        Geometry.cpp \
	Quaternion.cpp		\
	Dodecahedron.cpp	\
	Cylinder.cpp		\
	Boundary.cpp		\
	ComplexNumber.cpp	\
	FourierTransformable.cpp \
	FourierTransformable2D.cpp \

DISLOCATION_SOURCES = DislocationNode.cpp \
			DislocationSegment.cpp \
			SurfaceSegment.cpp

SOURCES = 	SurfaceLoop.cpp \
		SurfaceProfile.cpp \
		SurfaceDeformer.cpp 


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
        
        


