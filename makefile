SRCDIR = src
INCDIR = inc
OBJDIR = obj
BINDIR = bin

MAKEFILE_NAME = ${firstword ${MAKEFILE_LIST}}# makefile name

CXX=g++
CXXFLAGS=-I$(INCDIR) -std=c++14 -Wall -Wextra -MMD -g #-v

rm=rm -f

target1 = test
objnames1 = test.o Genome.o Parameters.o Random.o Phenotype.o TaskFunctor.o BaseFunctor.o InnovationTracker.o \
			FitnessFunctor.o LibNeat.o NEAT.o SimpleFitnessFunctor.o

# target2 = randomcode
# objnames2 = BFIInterpretter.o RandomCode.o

fulltarget1 = $(target1:%=$(BINDIR)/%)
# fulltarget2 = $(target2:%=$(BINDIR)/%)
targets = $(fulltarget1) # $(fulltarget2)

objects1 = $(objnames1:%=$(OBJDIR)/%)
# objects2 = $(objnames2:%=$(OBJDIR)/%)
objects = $(objects1) # $(objects2)

depends = ${objects:.o=.d}

.PHONY: clean all

all: $(targets)

$(fulltarget1) : $(objects1)
	$(CXX) $(CXXFLAGS) $^ -o $@

# $(fulltarget2) : $(objects2)
# 	$(CXX) $(CXXFLAGS) $^ -o $@

$(objects): $(OBJDIR)/%.o : $(SRCDIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

${objects} : ${MAKEFILE_NAME}

-include ${depends}

clean:
	$(rm) $(objects)
	@echo "Object files removed!"
	$(rm) $(depends)
	@echo "Depend files removed!"
	$(rm) $(targets)
	@echo "Executables removed!"
