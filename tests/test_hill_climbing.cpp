#include <array>
#include <cstdio>

#include "test_utils.h"
#include "sbox_properties.h"
#include "hill_climbing.h"
#include "utils.h"

template <typename T>
void log_data(
	sbgen::shared_info_t<T>& params,
	sbgen::hill_climbing_info_t<T>& info
	)
{
	std::cout<<"target s-box found!"<<std::endl;
	std::cout<<"iteration: "<<params.iteration.load()<<std::endl;
	return;
}

int test_hill_climbing_with_whs1() 
{
	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = true;
	info.use_log_function = true;
	info.default_log_output = false;
	info.log_good_nl = true;
	info.log_better_sbox = true;
	info.log_good_nl_function = log_data<double>;
	info.log_better_sbox_function = [](auto& params, auto& info) -> void
		{
			std::cout<<"better sbox found!"<<std::endl;
			std::cout<<"iteration: "<<params.iteration.load()<<std::endl;
			return;
		};

	info.try_per_thread = 1000000;
	info.max_frozen_count = 100000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);

	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::hill_climbing<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())!=102) {
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_hill_climbing_with_whs2() 
{

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 10;
	info.max_frozen_count = 100000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 106);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::hill_climbing<double>(info);

	if (sbox.has_value()) {
		return REJECT_TEST;
	}
	
	return ACCEPT_TEST;
}

int test_hill_climbing_with_whs3() 
{

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_count = 1000000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);
	setup_property( &info, SBGEN_DELTA_UNIFORMITY, 8);
	setup_property( &info, SBGEN_ALGEBRAIC_IMMUNITY, 3);
	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::hill_climbing<double>(info);

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


int test_hill_climbing_with_wcf1() 
{

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_count = 100000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::wcf<double>;
	//info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::hill_climbing<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())!=102) {
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_hill_climbing_with_wcf2() 
{

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 10;
	info.max_frozen_count = 100000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 106);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::wcf<double>;
	//info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::hill_climbing<double>(info);

	if (sbox.has_value()) {
		return REJECT_TEST;
	}
	
	return ACCEPT_TEST;
}

int test_hill_climbing_with_wcf3() 
{

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_count = 1000000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);
	setup_property( &info, SBGEN_DELTA_UNIFORMITY, 8);
	setup_property( &info, SBGEN_ALGEBRAIC_IMMUNITY, 3);
	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::wcf<double>;
	//info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	auto sbox = sbgen::hill_climbing<double>(info);

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

int test_hill_climbing_with_pcf1() 
{

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_count = 100000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::pcf<double>;
	info.cost_data.reset(new sbgen::pcf_function_data_t(5));

	auto sbox = sbgen::hill_climbing<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}
	
	if(sbgen::properties::nonlinearity(sbox.value())!=102) {
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_hill_climbing_with_pcf2() 
{

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 10;
	info.max_frozen_count = 100000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 106);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::pcf<double>;
	info.cost_data.reset(new sbgen::pcf_function_data_t(5));

	auto sbox = sbgen::hill_climbing<double>(info);

	if (sbox.has_value()) {
		return REJECT_TEST;
	}
	
	return ACCEPT_TEST;
}

int test_hill_climbing_with_pcf3() 
{

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = true;

	info.try_per_thread = 1000000;
	info.max_frozen_count = 1000000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);
	setup_property( &info, SBGEN_DELTA_UNIFORMITY, 8);
	setup_property( &info, SBGEN_ALGEBRAIC_IMMUNITY, 3);
	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::pcf<double>;
	info.cost_data.reset(new sbgen::pcf_function_data_t(5));

	auto sbox = sbgen::hill_climbing<double>(info);

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

int test_hill_climbing_with_cf1() 
{

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 8;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_count = 100000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 104);

	
	info.use_random_seed = true;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::cf1<double>;
	info.cost_data.reset(new sbgen::cf1_function_data_t(12, 32, 0));

	auto sbox = sbgen::hill_climbing<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}

	if(sbgen::properties::nonlinearity(sbox.value())!=104) {
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_hill_climbing_with_cf2_1() 
{

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_count = 100000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 104);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;

	info.cost_function = sbgen::cf2<double>;
	info.cost_data.reset(new sbgen::cf2_function_data_t(1, 32, 32));

	auto sbox = sbgen::hill_climbing<double>(info);

	if (!sbox.has_value()) {
		return REJECT_TEST;
	}

	if(sbgen::properties::nonlinearity(sbox.value())!=104) {
		return REJECT_TEST;
	}
	
	//PRINT_SBOX(sbox.value());
	return ACCEPT_TEST;
}

int test_hill_climbing_with_cf2_2() 
{

	sbgen::hill_climbing_info_t<uint64_t> info;

	info.thread_count = 1;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_count = 100000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 104);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;

	info.cost_function = sbgen::cf2<uint64_t>;
	info.cost_data.reset(new sbgen::cf2_function_data_t(1, 32, 32));

	auto sbox = sbgen::hill_climbing<uint64_t>(info);

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
	int total_test_count = 12;
	int accepted_test_count = 0;
	
	run_test(test_hill_climbing_with_whs1, 
		"test hill climbing with whs 1\n", accepted_test_count);
	run_test(test_hill_climbing_with_whs2,
		"test hill climbing with whs 2\n", accepted_test_count);
	run_test(test_hill_climbing_with_whs3, "test hill climbing with whs 3\n", accepted_test_count);
	
	run_test(test_hill_climbing_with_wcf1,
		"test hill climbing with wcf 1\n", accepted_test_count);
	run_test(test_hill_climbing_with_wcf2,
		"test hill climbing with wcf 2\n", accepted_test_count);
	run_test(test_hill_climbing_with_wcf3, "test hill climbing with wcf 3\n", accepted_test_count);
	
	run_test(test_hill_climbing_with_pcf1,
		"test hill climbing with pcf 1\n", accepted_test_count);
	run_test(test_hill_climbing_with_pcf2,
		"test hill climbing with pcf 2\n", accepted_test_count);
	run_test(test_hill_climbing_with_pcf3, "test hill climbing with pcf 3\n", accepted_test_count);
	run_test(test_hill_climbing_with_cf1,
		"test hill climbing with cf1\n", accepted_test_count);
	run_test(test_hill_climbing_with_cf2_1,
		"test hill climbing with cf2 1\n", accepted_test_count);
	run_test(test_hill_climbing_with_cf2_2,
		"test hill climbing with cf2 2\n", accepted_test_count);
	
	printf("%d/%d test passed\n", accepted_test_count, total_test_count);
	return !(total_test_count==accepted_test_count);
}
