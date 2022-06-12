// Simulated annealing implementation
#ifndef  _SBGEN_SIMULATED_ANNEALING_H_
#define _SBGEN_SIMULATED_ANNEALING_H_

#include <thread>
#include <functional>
#include <optional>
#include <random>
#include <vector>
#include <iostream>

#include "cost_function.h"
#include "generator.h"

namespace sbgen {

/**
 * @brief Simulated annealing method parameters
 **/
template <typename T>
struct simulated_annealing_info_t : public properties_info_t{
	// Thread count
	int32_t thread_count;
	// maximal try count in one thread
	int32_t try_per_thread;
	// maximal outer loops count
	int32_t max_outer_loops;
	// maximal inner loops count
	int32_t max_inner_loops;
	// maximal frozen loops count
	int32_t max_frozen_outer_loops;
	// visibility of computation process
	bool is_log_enabled;
    
	// annealing parameters
	double initial_temperature;
	double alpha_parameter;

	// cost function and cost function data
	std::unique_ptr< cost_function_data_t> cost_data;
	std::function< cost_info_t<T>(cost_function_data_t*, std::array<uint8_t, 256>)> cost_function;
}; // simulated_annealing_info_t


/**
  * @brief One-thread simulated annealing generator
  *
  * Simulated annealing for one thread
  *
  * @param params
  *    shared data
  * @param info
  *    simulated annealing parameters
  * @param id
  *    thread id (used for randomness)
  *@returns 
  *   s-box with target properties  in params
  */
template<typename T>
void simulated_annealing_thread_function (
	shared_info_t<T>& params, 
	simulated_annealing_info_t<T>& info, 
	int id
) 
{
	std::random_device rd;
	unsigned seed;
	bool      accept_in_this_loop = false;
	
	if (info.use_random_seed) 
	{
		seed = rd();
	} else 
	{
		seed = info.seed^id;
	}
		
	std::mt19937 gen(seed);
	std::uniform_int_distribution<int> distrib(0, 255);
	std::uniform_real_distribution<> dis(0.0, 1.0);
	std::array<uint8_t, 256> new_sbox;
	cost_info_t<T> new_cost;
	double    current_temperature = info.initial_temperature;

	for (int32_t x = 0;x < info.max_outer_loops;x++) {
		accept_in_this_loop = false;
		
		for (int32_t y = 0;y < info.max_inner_loops;y++) {
			params.iteration.fetch_add(1);
			
			int pos_1 = 0;
			int pos_2 = 0;

			params.sbox_mutex.lock();
			if (params.is_found) {
				params.sbox_mutex.unlock();
				return;
			}
			
			new_sbox = params.best_sbox;
			params.sbox_mutex.unlock();

			while (pos_1 == pos_2) {
				pos_1 = distrib(gen);
				pos_2 = distrib(gen);
			}
			std::swap(new_sbox[pos_1], new_sbox[pos_2]);

			new_cost = info.cost_function(info.cost_data.get(), new_sbox);
			bool is_good_nl = new_cost.nonlinearity >= info.target_properties[SBGEN_NONLINEARITY];
			while (is_good_nl) {
			
				if (info.is_log_enabled)
					std::cout << "cost=" << params.best_cost.cost << "	NL=" << params.best_cost.nonlinearity << "	iteration=" << params.iteration.load() << std::endl;
			
				if(!generator_utils::check_additional_properties(static_cast<properties_info_t>(info),new_sbox))
					break;
				params.sbox_mutex.lock();
				params.best_sbox = new_sbox;
				params.best_cost = new_cost;
				params.is_found = true;
				params.sbox_mutex.unlock();
				return;
			}
		
			params.sbox_mutex.lock();
			double cost_diff = (double)(new_cost.cost - params.best_cost .cost);

			if ((cost_diff <= 0 && new_cost.nonlinearity >= params.best_cost.nonlinearity) 
				|| new_cost.nonlinearity > params.best_cost.nonlinearity ) {
				params.best_sbox = new_sbox;
				params.best_cost = new_cost;
				accept_in_this_loop = true;

				if (info.is_log_enabled)
					std::cout << "cost=" << params.best_cost.cost << "	NL=" << params.best_cost.nonlinearity << "	iteration=" << params.iteration.load() << std::endl;
			} else {
				double u = dis(gen);
				if (u < exp(-cost_diff / current_temperature)) {
					params.best_sbox = new_sbox;
					params.best_cost = new_cost;
					accept_in_this_loop = true;
				}
			}
			
			if(params.frozen_count/info.thread_count >= info.max_frozen_outer_loops) {
				params.sbox_mutex.unlock();
				break;
			}
			
			if (accept_in_this_loop == false) {
				params.frozen_count++;
				if (params.frozen_count/info.thread_count >= info.max_frozen_outer_loops)
					break;
			} else
				params.frozen_count = 0;
			params.sbox_mutex.unlock();
			
			current_temperature *= info.alpha_parameter;

		}
		
	}
	return;
}


/**
  * @brief Mult-thread simulated annealing generator
  *
  * Simulated annealing generator
  *
  * @param info
  *    simulated annealing parameters
  *@returns 
  *   s-box with target properties 
  */
template<typename T>
std::optional<std::array<uint8_t, 256>> simulated_annealing (simulated_annealing_info_t<T>& info) {
	shared_info_t<T> thread_data;
	std::random_device rd;
	
	if (info.use_random_seed) 
	{
		info.seed = rd();
	}
	
	std::mt19937 gen(info.seed);
	std::uniform_int_distribution<int> distrib(0, 255);
	std::vector<std::thread> workers;


	for (int i = 0; i <= 255; i++) {
		thread_data.best_sbox[i] = i;
	}

	for (int i = 255; i > 0; i--) {
		int j = distrib(gen) % i;
		std::swap(thread_data.best_sbox[i], thread_data.best_sbox[j]);
	}

	thread_data.is_found = false;
	thread_data.frozen_count = 0;
	thread_data.iteration.store(0);
	thread_data.best_cost = info.cost_function(info.cost_data.get(), thread_data.best_sbox);

	for (int i = 0; i < info.thread_count; i++) {
		workers.push_back(
			std::thread(&simulated_annealing_thread_function<T>,
				std::ref(thread_data), std::ref(info), i)
		);
	}

	for (std::thread& worker : workers)
		worker.join();
	workers.clear();

	if (thread_data.is_found) {
		return thread_data.best_sbox;
	}

	return {};
}
	
}; // namespace sbgen

#endif // _SBGEN_SIMULATED_ANNEALING_H_
