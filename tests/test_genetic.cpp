#include <array>
#include <cstdio>

#include "test_utils.h"
#include "sbox_properties.h"
#include "genetic.h"
#include "genetic_internals.h"
#include "utils.h"

int test_hill_climbing_with_whs1() 
{

	sbgen::genetic_info_t<double> info;

	info.thread_count = 8;
	info.is_log_enabled = true;

	info.mutants_per_parent = 10;
	info.selection_count = 10;
	info.iterations_count = 15000;
	info.initial_population_count = 100;
	
	info.crossover_count = 0;
	info.child_per_parent = 0;
	info.crossover_method = 
	info.use_crossover = false;
	
	setup_property( &info, SBGEN_NONLINEARITY, 104);

	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));
	
	info.selection_method = sbgen::selectors::basic_selection<double>;
	info.comparator.comparator = sbgen::comparators::less<double>;

	auto sbox = sbgen::genetic<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())!=104) {
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int main()
{
	int total_test_count = 1;
	int accepted_test_count = 0;
	
	run_test(test_hill_climbing_with_whs1, "test genetic with whs 1\n", accepted_test_count);

	printf("%d/%d test passed\n", accepted_test_count, total_test_count);
	return !(total_test_count==accepted_test_count);
}
