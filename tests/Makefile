CXXFLAGS += -std=c++17 -O3 -Wall -pthread
CPPFLAGS += -I ../include

TESTS = test_properties  test_hill_climbing test_simulated_annealing

all:  $(TESTS)


test_properties:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test_properties.cpp  -o $@
test_hill_climbing:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test_hill_climbing.cpp  -o $@
test_simulated_annealing:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test_simulated_annealing.cpp  -o $@
	
clean:
	rm $(TESTS)
run:
	@echo "test properties"
	@./test_properties
	@echo "test hill climbing"
	@./test_hill_climbing
	@echo "test simulated annealing"
	@./test_simulated_annealing

