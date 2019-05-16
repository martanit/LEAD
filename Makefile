#
# LEDy, Loop-extrusion molecular dynamics. Code is hosted on GitHub.
#
# File: Makefile
#
# Authors:
#  Martina Crippa 				<martina.crippa2@studenti.unimi.it>
#
################################################################################

# configure multithreaded compilation
OS := $(shell uname)
ifeq ($(OS), Linux)
	export MAKEFLAGS="-j $(grep -c ^processor /proc/cpuinfo)"
else ifeq ($(OS), Darwin)
	export MAKEFLAGS="-j $(sysctl-n hw.ncpu)"
endif

# search in VPATH for source file
VPATH=./source/
# place object file in OBJPATH
OBJPATH=./obj/

CXX := $(GCC47_BINDIR)$(CXX)

TARGET := LEDy
OBJS := $(patsubst %.o,$(OBJPATH)%.o, parameters.o polymer.o potential.o integrator.o dynamics.o )

DEBUG := -g
WARNING := -Wall -Wextra

CXXFLAGS := $(CXXFLAGS) -std=c++17 -O2
LDFLAGS := -lpthread

all: $(TARGET)

LEDy: $(OBJS) $(OBJPATH)main.o
	$(CXX) -o $@ $(OBJS) $(OBJPATH)main.o $(CXXFLAGS) $(LDFLAGS)

$(OBJPATH)main.o: main.cpp $(OBJS)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJPATH)%.o: %.cpp %.h
ifeq ($(wildcard $(OBJPATH)*),) #search for obj path, create it if it doesn't exist
	@mkdir -p $(OBJPATH)
endif
	$(CXX) -c $< -o $@ $(CXXFLAGS)

.PHONY: run clean

run:
	./$(TARGET)

clean:
	/bin/rm -f $(OBJPATH)*.o
	/bin/rm -f $(TARGET)

# vim: set noexpandtab:

