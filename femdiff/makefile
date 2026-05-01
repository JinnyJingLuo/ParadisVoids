# Declaration of variables
CC = mpicxx
EZMATH_DIR = ../ezmath/
FLAGS = -I$(EZMATH_DIR) -DPARALLEL -vec_report0
OPTIMIZATION = -O3
EXEC = femdiff

EZMATH_SOURCES = Point.cpp				\
	Tools.cpp				\
	Vector.cpp				\
	Matrix.cpp				\
	Randomizer.cpp

SOURCES =  FEMDiff.cpp

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



