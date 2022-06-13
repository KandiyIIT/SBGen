// Hill climbing implementation
#ifndef _HILL_CLIMBING_H_
#define _HILL_CLIMBING_H_

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
 * @brief Hill climbing method parameters
 **/
template <typename T>
struct hill_climbing_info_t : public properties_info_t{
    // thread count
	int32_t thread_count;
    // maximal try count in one thread
	int32_t try_per_thread;
    // maximal frozen loops count
	int32_t max_frozen_count;
    // visibility of computation process
	bool is_log_enabled;

    // cost function and cost function data
	std::unique_ptr< cost_function_data_t> cost_data;
	std::function< cost_info_t<T>(cost_function_data_t*, std::array<uint8_t, 256>)> cost_function;
}; // struct hill_climbing_info_t 


/**
  * @brief One-thread hiil climbing generator
  *
  * Hill Climbing for one thread
  *
  * @param params
  *    shared data
  * @param info
  *    hill climbing parameters
  * @param id
  *    thread id (used for randomness)
  *@returns 
  *   s-box with target properties  in params
  */
template<typename T>
void hill_climbing_thread_function (
	shared_info_t<T>& params, 
	hill_climbing_info_t<T>& info, int id
) 
{
	std::random_device rd;
	unsigned seed;
	
	if (info.use_random_seed) 
	{
		seed = rd();
	} else 
	{
		seed = info.seed^id;
	}
		
	std::mt19937 gen(seed);
	std::uniform_int_distribution<int> distrib(0, 255);
	std::array<uint8_t, 256> new_sbox;
	cost_info_t<T> new_cost;

	for (int i = 0; i < info.try_per_thread; i++) {
		params.iteration.fetch_add(1);

		int pos_1 = 0;
		int pos_2 = 0;

		params.sbox_mutex.lock();
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
		
		if(!is_good_nl){
			params.sbox_mutex.lock();
			params.frozen_count++;
			if (params.frozen_count > info.max_frozen_count) {
				if (info.is_log_enabled)
					std::cout << "iteration=" << params.iteration.load() << "	Search stopped: frozen_count > max_frozen_count" << std::endl;
				params.sbox_mutex.unlock();
				return;
			}
			params.sbox_mutex.unlock();
		}

		params.sbox_mutex.lock();
		if (params.is_found) {
			params.sbox_mutex.unlock();
			return;
		}

		if ((new_cost.cost <= params.best_cost.cost && new_cost.nonlinearity >= params.best_cost.nonlinearity) 
			|| new_cost.nonlinearity > params.best_cost.nonlinearity ) {
			params.best_sbox = new_sbox;
			params.best_cost = new_cost;
			params.frozen_count = 0;

			if (info.is_log_enabled)
				std::cout << "cost=" << params.best_cost.cost << "	NL=" << params.best_cost.nonlinearity << "	iteration=" << params.iteration.load() << std::endl;

		}
		params.sbox_mutex.unlock();
	}
	return;
}

/**
  * @brief Mult-thread hiil climbing generator
  *
  * Hill Climbing generator
  *
  * @param info
  *    hill climbing parameters
  *@returns 
  *   s-box with target properties 
  */
template<typename T>
std::optional<std::array<uint8_t, 256>> hill_climbing(hill_climbing_info_t<T>& info) {
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
			std::thread(&hill_climbing_thread_function<T>,
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

}; // _HILL_CLIMBING_H_

#endif // _HILL_CLIMBING_H_
