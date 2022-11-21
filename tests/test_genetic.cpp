#include <array>
#include <cstdio>

#include "test_utils.h"
#include "sbox_properties.h"
#include "genetic.h"
#include "genetic_internals.h"
#include "utils.h"

int test_genetic_with_basic_selection() 
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
	//info.crossover_method = 
	info.use_crossover = false;
	
	setup_property( &info, SBGEN_NONLINEARITY, 104);
	setup_property( &info, SBGEN_DELTA_UNIFORMITY, 8);
	setup_property( &info, SBGEN_ALGEBRAIC_IMMUNITY, 3);

	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));
	
	info.selection_method = sbgen::selectors::basic_selection<double>;
	info.comparator.comparator = sbgen::comparators::less_nl<double>;

	auto sbox = sbgen::genetic<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())!=104) {
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	if(sbgen::properties::nonlinearity(sbox.value())<104) {
		printf("Nl=%d\n",sbgen::properties::nonlinearity(sbox.value()));
		return REJECT_TEST;
	}
	
	if(sbgen::properties::delta_uniformity(sbox.value())>8) {
		printf("DU=%d\n",sbgen::properties::delta_uniformity(sbox.value()));
		return REJECT_TEST;
	}
	
	if(sbgen::properties::algebraic_immunity(sbox.value())<3) {
		printf("AI=%d\n",sbgen::properties::algebraic_immunity(sbox.value()));
		return REJECT_TEST;
	}
	return ACCEPT_TEST;
}

int test_genetic_with_rank_selection() 
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
	//info.crossover_method = 
	info.use_crossover = false;
	
	setup_property( &info, SBGEN_NONLINEARITY, 104);

	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));
	
	info.selection_method = sbgen::selectors::rank_sequential_selection<double>;
	info.comparator.comparator = sbgen::comparators::less_nl<double>;

	auto sbox = sbgen::genetic<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())<104) {
		printf("Nl=%d\n",sbgen::properties::nonlinearity(sbox.value()));
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_genetic_with_roulette_wheel_selection() 
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
	//info.crossover_method = 
	info.use_crossover = false;
	
	setup_property( &info, SBGEN_NONLINEARITY, 104);

	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));
	
	info.selection_method =
		sbgen::selectors::roulette_wheel_sequential_selection<double>;
	info.comparator.comparator = sbgen::comparators::less_nl<double>;

	auto sbox = sbgen::genetic<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())<104) {
		printf("Nl=%d\n",sbgen::properties::nonlinearity(sbox.value()));
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_genetic_with_pmx_crossover() 
{

	sbgen::genetic_info_t<double> info;

	info.thread_count = 8;
	info.is_log_enabled = true;

	info.mutants_per_parent = 10;
	info.selection_count = 10;
	info.iterations_count = 15000;
	info.initial_population_count = 100;
	
	info.crossover_count = 50;
	info.child_per_parent = 1;
	info.crossover_method = sbgen::cossovers::pmx;
	info.use_crossover = true;
	
	setup_property( &info, SBGEN_NONLINEARITY, 104);

	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));
	
	info.selection_method =
		sbgen::selectors::basic_selection<double>;
		//roulette_wheel_sequential_selection<double>;
	info.comparator.comparator = sbgen::comparators::less_nl<double>;

	auto sbox = sbgen::genetic<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())<104) {
		printf("Nl=%d\n",sbgen::properties::nonlinearity(sbox.value()));
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_genetic_with_cyclic_crossover() 
{

	sbgen::genetic_info_t<double> info;

	info.thread_count = 8;
	info.is_log_enabled = true;

	info.mutants_per_parent = 10;
	info.selection_count = 10;
	info.iterations_count = 15000;
	info.initial_population_count = 100;
	
	info.crossover_count = 50;
	info.child_per_parent = 1;
	info.crossover_method = sbgen::cossovers::cycle;
	info.use_crossover = true;
	
	setup_property( &info, SBGEN_NONLINEARITY, 104);

	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));
	
	info.selection_method =
		sbgen::selectors::basic_selection<double>;
		//roulette_wheel_sequential_selection<double>;
	info.comparator.comparator = sbgen::comparators::less_nl<double>;

	auto sbox = sbgen::genetic<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())<104) {
		printf("Nl=%d\n",sbgen::properties::nonlinearity(sbox.value()));
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_genetic_106() 
{

	sbgen::genetic_info_t<double> info;

	info.thread_count = 8;
	info.is_log_enabled = true;

	info.mutants_per_parent = 10;
	info.selection_count = 10;
	info.iterations_count = 15000;
	info.initial_population_count = 100;
	
	info.crossover_count = 50;
	info.child_per_parent = 1;
	info.crossover_method = sbgen::cossovers::cycle;
	info.use_crossover = true;
	
	setup_property( &info, SBGEN_NONLINEARITY, 106);

	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));
	
	info.selection_method =
		sbgen::selectors::roulette_wheel_sequential_selection<double>;
	info.comparator.comparator = sbgen::comparators::less_nl<double>;

	auto sbox = sbgen::genetic<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())<106) {
		printf("Nl=%d\n",sbgen::properties::nonlinearity(sbox.value()));
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int main()
{
	int total_test_count = 5;
	int accepted_test_count = 0;
	//test_genetic_106();
	
	run_test(test_genetic_with_basic_selection, 
		"test genetic with basic selection\n", accepted_test_count);
	run_test(test_genetic_with_rank_selection, 
		"test genetic with rank selection\n", accepted_test_count);
	run_test(test_genetic_with_roulette_wheel_selection, 
		"test genetic with roulette wheel selection\n", accepted_test_count);
	run_test(test_genetic_with_pmx_crossover,
		"test genetic with pmx crossover\n", accepted_test_count);
	run_test(test_genetic_with_cyclic_crossover,
		"test genetic with cyclic crossover\n", accepted_test_count);

	printf("%d/%d test passed\n", accepted_test_count, total_test_count);
	return !(total_test_count==accepted_test_count);
}
