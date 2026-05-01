# Declaration of variables
CC = mpicxx
EZMATH_DIR = ../ezmath/
FLAGS = -I$(EZMATH_DIR) #-DPARALLEL
OPTIMIZATION = -g -O3
EXEC = casfem

EZMATH_SOURCES = AxisAlignedBoundingBox.cpp	\
	MathServices.cpp			\
	Matrix.cpp				\
	CartesianOrthogonalCoordinateSystem.cpp	\
	Plane.cpp				\
	Point.cpp				\
	PointPopulatedAxisAlignedBoundingBox.cpp\
	QuadPatch.cpp				\
	Randomizer.cpp				\
	Segment.cpp				\
	GenericNode.cpp				\
	Tools.cpp				\
	Patch.cpp				\
	TriPatch.cpp				\
	Line.cpp				\
	Vector.cpp				\
	SparseMatrix.cpp		\
	GeometricComponent.cpp		\
	Geometry.cpp \
	Block.cpp	\
	Cylinder.cpp	\

SOURCES = FEMMesh.cpp \
		FEMDegreeOfFreedom.cpp \
		FEMPlasticityState.cpp \
		FEMNode.cpp \
		FEMSolidNode.cpp \
		FEMPotentialNode.cpp \
		FEMThermoMechanicalNode.cpp \
		FEMLoad.cpp \
		FEMConstantLoad.cpp \
		FEMTimeSinusoidalLoad.cpp \
		FEMTimeLinearLoad.cpp \
		FEMXTorsionLoad.cpp \
		FEMYTorsionLoad.cpp \
		FEMAxialTimeLinearTotalStrainLoad.cpp \
		FEMAxialTimeSinusoidalTotalStrainLoad.cpp \
		FEMStepLoad.cpp \
		FEMPulseTrainLoad.cpp \
		FEMMaterial.cpp \
		FEMLinearElasticIsotropicMaterial.cpp \
		FEMJ2PlasticIsotropicMaterial.cpp \
		FEMGaussPoint.cpp \
		FEMElementGeometry.cpp \
		FEMHexahedralElement.cpp \
		FEMElement.cpp \
		FEMSolidElement.cpp \
		FEMPotentialElement.cpp \
		FEMThermoMechanicalElement.cpp \
		FEMSolver.cpp \
		FEMSolidSolver.cpp \
		FEMThermoMechanicalSolver.cpp \
		FEMImplicitDynamicsSolidSolver.cpp \
		FEMPotentialSolver.cpp \
		FEMExplicitDynamicsSolidSolver.cpp \
		FEMImplicitDynamicsThermoMechanicalSolver.cpp \
		FEMBoundaryElementFace.cpp \
		MainDataStructure.cpp \
		CASFEM.cpp

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



