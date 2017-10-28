SRCDIR = src/
INCDIR = ./include
OBJDIR = objs/
SIMDIR = sims/
MATHGEODIR = /Users/dhoff/Programs/MathGeoLib/MathGeoLib/src/

OBJECTS := objs/material.o objs/planeSegment.o objs/plane.o objs/detector.o objs/plf.o objs/frag.o objs/tele.o objs/array.o objs/sle.o objs/decay.o objs/loss.o objs/mScat.o objs/fragment.o objs/random.o objs/tower.o objs/frame.o objs/MomDist.o objs/moscow.o objs/ring.o objs/ttt.o objs/plane_det.o objs/detectorV2.o

#ALLOBJECTS = $(wildcard $(SRCDIR)*.cpp)

ALLOBJECTS := $(patsubst $(SRCDIR)%.cpp,$(OBJDIR)%.o,$(wildcard $(SRCDIR)*.cpp))

CFLAGS= -c -g -std=c++11 -I$(shell root-config --incdir) -I$(INCDIR) -I$(MATHGEODIR)
COMPILER= clang++
LINKOPTION = $(shell root-config --libs)  -L$(MATHGEODIR)../osx_make/ -lMathGeoLib

all: simCa34 simK33 simO26

simO26: $(OBJDIR)simO26.o $(OBJECTS)
	$(COMPILER) -o sims/simO26 $(OBJDIR)simO26.o $(OBJECTS) $(LINKOPTION)

simK33: $(OBJDIR)simK33.o $(OBJECTS)
	$(COMPILER) -o sims/simK33 $(OBJDIR)simK33.o $(OBJECTS) $(LINKOPTION)

simCa34: $(OBJDIR)simCa34.o $(OBJECTS)
	$(COMPILER) -o sims/simCa34 $(OBJDIR)simCa34.o $(OBJECTS) $(LINKOPTION)

$(ALLOBJECTS) : $(OBJDIR)%.o : $(SRCDIR)%.cpp
	$(COMPILER) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJDIR)*.o

