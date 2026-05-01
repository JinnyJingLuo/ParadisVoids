# Declaration of variables
CC = mpicxx
EZMATH_DIR = ../ezmath/
FLAGS = -vec_report0 -I$(EZMATH_DIR)
OPTIMIZATION = -O3
EXEC = newDD

EZMATH_SOURCES = AxisAlignedBoundingBox.cpp	\
	MathServices.cpp			\
	Matrix.cpp				\
	CartesianOrthogonalCoordinateSystem.cpp	\
	CategorizedGenericNode.cpp		\
	Patch.cpp				\
	Plane.cpp				\
	Point.cpp				\
	PointPopulatedAxisAlignedBoundingBox.cpp\
	QuadPatch.cpp				\
	Randomizer.cpp				\
	Segment.cpp				\
	SparseMatrix.cpp			\
	GenericNode.cpp				\
	Tools.cpp				\
	TriPatch.cpp				\
	Line.cpp				\
	Geometry.cpp			\
	Vector.cpp

SOURCES = DislocationNodeTag.cpp			\
	DislocationNodeTagMap.cpp		\
	DislocationNode.cpp			\
	DislocationSegment.cpp			\
	MainDataStructure.cpp			\
	ChangeConnectionOperation.cpp		\
	TopologicalOperation.cpp		\
	TopologyServer.cpp			\
	DynamicsServer.cpp			\
	SurfaceServer.cpp			\
	BalancingServer.cpp			\
	DislocationSimulator.cpp		

OBJECTS = $(SOURCES:.cpp=.o)

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
        
