#include <array>
#include <cstdio>

#include "test_utils.h"
#include "sbox_properties.h"
#include "simulated_annealing.h"
#include "utils.h"



int test_simulated_annealing_with_whs1() 
{

	sbgen::simulated_annealing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_outer_loops = 100000;
    info.max_outer_loops = 10000;
    info.max_inner_loops = 1000;
    info.initial_temperature = 100;
    info.alpha_parameter = 0.99;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::simulated_annealing<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())!=102) {
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_simulated_annealing_with_whs2() 
{

	sbgen::simulated_annealing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 10;
	info.max_frozen_outer_loops = 100000;
    info.max_outer_loops = 10;
    info.max_inner_loops = 1;
    info.initial_temperature = 100;
    info.alpha_parameter = 0.99;
	
	setup_property( &info, SBGEN_NONLINEARITY, 106);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::simulated_annealing<double>(info);

	if (sbox.has_value()) {
		return REJECT_TEST;
	}
	
	return ACCEPT_TEST;
}

int test_simulated_annealing_with_whs3() 
{

	sbgen::simulated_annealing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_outer_loops = 100000;
    info.max_outer_loops = 10000;
    info.max_inner_loops = 1000;
    info.initial_temperature = 100;
    info.alpha_parameter = 0.99;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);
	setup_property( &info, SBGEN_DELTA_UNIFORMITY, 8);
	setup_property( &info, SBGEN_ALGEBRAIC_IMMUNITY, 3);
	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::simulated_annealing<double>(info);

	if (!sbox.has_value()) {
		printf("not found!\n");
		return REJECT_TEST;
	}
	
	
	if(sbgen::properties::nonlinearity(sbox.value())<102) {
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


int test_simulated_annealing_with_wcf1() 
{

	sbgen::simulated_annealing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_outer_loops = 100000;
    info.max_outer_loops = 10000;
    info.max_inner_loops = 1000;
    info.initial_temperature = 100;
    info.alpha_parameter = 0.99;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::wcf<double>;
	//info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::simulated_annealing<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())!=102) {
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_simulated_annealing_with_wcf2() 
{

	sbgen::simulated_annealing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 10;
	info.max_frozen_outer_loops = 100000;
    info.max_outer_loops = 1;
    info.max_inner_loops = 1;
    info.initial_temperature = 100;
    info.alpha_parameter = 0.99;
	
	setup_property( &info, SBGEN_NONLINEARITY, 106);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::wcf<double>;
	//info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::simulated_annealing<double>(info);

	if (sbox.has_value()) {
		return REJECT_TEST;
	}
	
	return ACCEPT_TEST;
}

int test_simulated_annealing_with_wcf3() 
{

	sbgen::simulated_annealing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_outer_loops = 100000;
    info.max_outer_loops = 10000;
    info.max_inner_loops = 1000;
    info.initial_temperature = 100;
    info.alpha_parameter = 0.99;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);
	setup_property( &info, SBGEN_DELTA_UNIFORMITY, 8);
	setup_property( &info, SBGEN_ALGEBRAIC_IMMUNITY, 3);
	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::wcf<double>;
	//info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::simulated_annealing<double>(info);

	if (!sbox.has_value()) {
		printf("not found!\n");
		return REJECT_TEST;
	}
	
	
	if(sbgen::properties::nonlinearity(sbox.value())<102) {
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

int test_simulated_annealing_with_pcf1() 
{

	sbgen::simulated_annealing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = true;

	info.try_per_thread = 1000000;
	info.max_frozen_outer_loops = 100000;
    info.max_outer_loops = 10000;
    info.max_inner_loops = 1000;
    info.initial_temperature = 1000000;
    info.alpha_parameter = 1;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::pcf<double>;
	info.cost_data.reset(new sbgen::pcf_function_data_t(7));

	auto sbox = sbgen::simulated_annealing<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())!=102) {
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_simulated_annealing_with_pcf2() 
{

	sbgen::simulated_annealing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 10;
	info.max_frozen_outer_loops = 100000;
    info.max_outer_loops = 1;
    info.max_inner_loops = 1;
    info.initial_temperature = 1000000;
    info.alpha_parameter = 1;
	
	setup_property( &info, SBGEN_NONLINEARITY, 106);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::pcf<double>;
	info.cost_data.reset(new sbgen::pcf_function_data_t(7));

	auto sbox = sbgen::simulated_annealing<double>(info);

	if (sbox.has_value()) {
		return REJECT_TEST;
	}
	
	return ACCEPT_TEST;
}

int test_simulated_annealing_with_pcf3() 
{

	sbgen::simulated_annealing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_outer_loops = 100000;
    info.max_outer_loops = 10000;
    info.max_inner_loops = 1000;
    info.initial_temperature = 1000000;
    info.alpha_parameter = 1;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);
	setup_property( &info, SBGEN_DELTA_UNIFORMITY, 8);
	setup_property( &info, SBGEN_ALGEBRAIC_IMMUNITY, 3);
	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::pcf<double>;
	info.cost_data.reset(new sbgen::pcf_function_data_t(7));

	auto sbox = sbgen::simulated_annealing<double>(info);

	if (!sbox.has_value()) {
		printf("not found!\n");
		return REJECT_TEST;
	}
	
	
	if(sbgen::properties::nonlinearity(sbox.value())<102) {
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

int main()
{
	int total_test_count = 4;
	int accepted_test_count = 0;
	
	run_test(test_simulated_annealing_with_whs1, "test simulated annealing with whs 1\n", accepted_test_count);
	run_test(test_simulated_annealing_with_whs2, "test simulated annealing with whs 2\n", accepted_test_count);
    //run_test(test_simulated_annealing_with_whs3, "test simulated annealing with whs 3\n", accepted_test_count);
	
	run_test(test_simulated_annealing_with_wcf1, "test simulated annealing with wcf 1\n", accepted_test_count);
	run_test(test_simulated_annealing_with_wcf2, "test simulated annealing with wcf 2\n", accepted_test_count);
	//run_test(test_simulated_annealing_with_wcf3, "test simulated annealing with wcf 3\n", accepted_test_count);
	
	//run_test(test_simulated_annealing_with_pcf1, "test simulated annealing with pcf 1\n", accepted_test_count);
	//run_test(test_simulated_annealing_with_pcf2, "test simulated annealing with pcf 2\n", accepted_test_count);
	//run_test(test_simulated_annealing_with_pcf3, "test simulated annealing with pcf 3\n", accepted_test_count);
	
	printf("%d/%d test passed\n", accepted_test_count, total_test_count);
	return !(total_test_count==accepted_test_count);
}
