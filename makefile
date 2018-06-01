CCX=g++
CXXFLAGS= -g -Wall -std=c++11 -w

DEP_FLAGS=-MMD -MP

SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
DEP=$(SRC:.cpp=.d)

CXXFLAGS+=$(DEP_FLAGS)

TARG=Project6

all: $(TARG) 

$(TARG): $(OBJ) 
	$(CCX) $^ $(LFLAGS) -o $@


.PHONY: clean

clean:
	rm $(OBJ) $(TARG) $(DEP)

-include $(DEP)