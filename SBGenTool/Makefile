CXXFLAGS += -std=c++17 -O3 -Wall -pthread
CPPFLAGS += -I ../include

TARGETS = sbgen

all:  $(TARGETS)

sbgen:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) SBGen.cpp  -o $@
	
clean:
	rm $(TARGETS)
