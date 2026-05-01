# Declaration of variables
CC = mpicxx
EZMATH_DIR = ../ezmath/
CASFEM_DIR = ../CASFEM/
DISLOCATION_DIR = dislocation_base/
FLAGS = -std=c++11 -I$(EZMATH_DIR) -I$(DISLOCATION_DIR) -I$(CASFEM_DIR) -DPARALLEL
OPTIMIZATION = -O2
EXEC = paradisJHU


EZMATH_SOURCES = AxisAlignedBoundingBox.cpp	\
	MathServices.cpp			\
	Matrix.cpp				\
	CartesianOrthogonalCoordinateSystem.cpp	\
	Patch.cpp				\
	Plane.cpp				\
	Point.cpp				\
	PointPopulatedAxisAlignedBoundingBox.cpp\
	QuadPatch.cpp				\
	Randomizer.cpp				\
	Segment.cpp				\
	GenericNode.cpp				\
	Tools.cpp				\
	TriPatch.cpp				\
	Line.cpp				\
	Vector.cpp				\
	Dodecahedron.cpp			\
	SparseMatrix.cpp			\
	Geometry.cpp				\
	Block.cpp				\
	Cylinder.cpp				\
	GeometricComponent.cpp			\
	Edge.cpp				\
	Face.cpp				\
	PolyhedronFace.cpp			\
	Polyhedron.cpp				\
	Curve.cpp

CASFEM_SOURCES = FEMMesh.cpp \
		FEMDegreeOfFreedom.cpp \
                FEMNode.cpp \
                FEMPotentialNode.cpp \
		FEMThermoMechanicalNode.cpp \
		FEMSolidNode.cpp \
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
                FEMElementGeometry.cpp \
                FEMHexahedralElement.cpp \
                FEMElement.cpp \
                FEMSolidElement.cpp \
		FEMPotentialElement.cpp \
		FEMThermoMechanicalElement.cpp \
		FEMThermoMechanicalSolver.cpp \
                FEMSolver.cpp \
                FEMSolidSolver.cpp \
                FEMPotentialSolver.cpp \
		FEMImplicitDynamicsSolidSolver.cpp \
                FEMExplicitDynamicsSolidSolver.cpp \
                FEMImplicitDynamicsThermoMechanicalSolver.cpp \
                FEMBoundaryElementFace.cpp \
		FEMGaussPoint.cpp \
		FEMPlasticityState.cpp \
		FEMJ2PlasticIsotropicMaterial.cpp \
		MainDataStructure.cpp \

DISLOCATION_SOURCES = DislocationNode.cpp \
			DislocationSegment.cpp \
			DislocationChain.cpp \
			SurfaceSegment.cpp \
			PrecipitateStructure.cpp

SOURCES = CellCharge.cpp  \
			Collision.cpp              \
				CommSendGhosts.cpp         \
				CommSendRemesh.cpp         \
				CommSendSecondaryGhosts.cpp \
				CommSendSegments.cpp       \
				CommSendVelocity.cpp       \
				Decomp.cpp                 \
				DeltaPlasticStrain.cpp     \
				DeltaPlasticStrain_FCC.cpp \
				DisableUnneededParams.cpp  \
	DLBfreeOld.cpp             \
				deWitInteraction.cpp       \
				DSMPI.cpp			\
	FindPreciseGlidePlane.cpp  \
				FixRemesh.cpp              \
				ForwardEulerIntegrator.cpp \
				FreeInitArrays.cpp         \
				GenerateOutput.cpp         \
				GetDensityDelta.cpp        \
				GetNewNativeNode.cpp       \
				GetNewGhostNode.cpp        \
				Heap.cpp                   \
				InitCellDomains.cpp        \
				InitCellNatives.cpp        \
				InitCellNeighbors.cpp      \
				InitRemoteDomains.cpp      \
      InitSendDomains.cpp        \
      Initialize.cpp             \
      InputSanity.cpp            \
      LoadCurve.cpp              \
      LocalSegForces.cpp         \
      ParadisMatrix.cpp                 \
      Meminfo.cpp                \
      Migrate.cpp                \
      MobilityLaw_FCC_0.cpp      \
      NodeForce.cpp              \
      NodeVelocity.cpp           \
      OsmoticForce.cpp           \
      ParadisFinish.cpp          \
      ParadisStep.cpp            \
      ParadisCrossSlipServer.cpp	\
	ParadisSurface.cpp	 \
	ParadisCylinderSurface.cpp \
	ParadisBlockSurface.cpp \
	ParadisTriangulatedSurface.cpp \
	ParadisGrainSurfaceServer.cpp	\
	ParadisExternalLoadServer.cpp \
	APBEventCommunicator.cpp \
	ParadisPrecipitateServer.cpp \
	ParadisCollisionServer.cpp \
	Param.cpp                  \
      Parse.cpp                  \
      PickScrewGlidePlane.cpp    \
      QueueOps.cpp               \
      ReadRestart.cpp            \
      RSDecomp.cpp               \
      RemapInitialTags.cpp       \
      Remesh.cpp                 \
      RemeshRule_2.cpp           \
      RemoveNode.cpp             \
      SortNativeNodes.cpp        \
      SortNodesForCollision.cpp  \
	Tecplot.cpp                \
      Topology.cpp               \
      Util.cpp                   \
      WriteDensityField.cpp      \
      WriteProp.cpp              \
      WriteRestart.cpp           \
	LoopFilter.cpp		\
	TwinPlaneCrossSlip.cpp \
      Main.cpp

ABSOLUTE_EZMATH_SOURCES = $(foreach source,$(EZMATH_SOURCES),$(EZMATH_DIR)$(source))
EZMATH_OBJECTS = $(ABSOLUTE_EZMATH_SOURCES:.cpp=.o)
ABSOLUTE_DISLOCATION_SOURCES = $(foreach source,$(DISLOCATION_SOURCES),$(DISLOCATION_DIR)$(source))
DISLOCATION_OBJECTS = $(ABSOLUTE_DISLOCATION_SOURCES:.cpp=.o)
ABSOLUTE_CASFEM_SOURCES = $(foreach source,$(CASFEM_SOURCES),$(CASFEM_DIR)$(source))
CASFEM_OBJECTS = $(ABSOLUTE_CASFEM_SOURCES:.cpp=.o)
OBJECTS = $(SOURCES:.cpp=.o)


# Main target
$(EXEC): $(EZMATH_OBJECTS) $(DISLOCATION_OBJECTS) $(CASFEM_OBJECTS) $(OBJECTS)
	$(CC) $(EZMATH_OBJECTS) $(DISLOCATION_OBJECTS) $(CASFEM_OBJECTS) $(OBJECTS) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(FLAGS) $(OPTIMIZATION) $< -o $@
#	$(CC) -c $(FLAGS) -O0 $< -o $@ -g

lint:
	clang-format --style="{SortIncludes: false}" -i *cpp *h dislocation_base/*cpp dislocation_base/*h

# To remove generated files
clean:
	rm -f $(EXEC) ./*.o
