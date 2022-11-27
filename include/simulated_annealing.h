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

namespace sbgen
{
	/**
	* @brief Simulated annealing method parameters
	**/
	template <typename T>
	struct simulated_annealing_info_t : public properties_info_t
	{
		using log_function_t =
			std::function<void(shared_info_t<T>&, simulated_annealing_info_t<T>&)>;
		int32_t thread_count;
		int32_t try_per_thread;
		int32_t max_outer_loops;
		int32_t max_inner_loops;
		int32_t max_frozen_outer_loops;

		bool is_log_enabled;
		bool use_log_function = false;
		bool default_log_output = true;
		bool log_good_nl = false;
		bool log_better_sbox = false;
    
		// annealing parameters
		double initial_temperature;
		double alpha_parameter;

		// functional parameters
		std::unique_ptr<cost_function_data_t> cost_data;
		cost_function_t<T> cost_function;
		comparator_t<T> comparator = comparators::less_nl<T>;
		log_function_t log_good_nl_function;
		log_function_t log_better_sbox_function;
	}; // simulated_annealing_info_t

	/**
	* @brief Simulated annealing log function type
	**/
	template <typename T>
	using simulated_annealing_log_function_t =
		typename simulated_annealing_info_t<T>::log_function_t;

	/**
	* @brief One-thread simulated annealing generator
	*
	* Simulated annealing for one thread
	*
	* @param params
	*    shared data
	* @param info
	*    simulated annealing parameters
	*@returns 
	*   s-box with target properties  in params
	*/
	template<typename T>
	void simulated_annealing_thread_function(shared_info_t<T>& params, 
		simulated_annealing_info_t<T>& info) 
	{
		bool accept_in_this_loop = false;
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> distrib(0, 255);
		std::uniform_real_distribution<> dis(0.0, 1.0);

		sbox_info<T> sbox;
		double    current_temperature = info.initial_temperature;

		for (int32_t x = 0;x < info.max_outer_loops;x++)
		{
			accept_in_this_loop = false;

			for (int32_t y = 0;y < info.max_inner_loops;y++)
			{
				params.iteration.fetch_add(1);

				int32_t pos_1 = 0;
				int32_t pos_2 = 0;
				while (pos_1 == pos_2)
				{
					pos_1 = distrib(gen);
					pos_2 = distrib(gen);
				}

				params.sbox_mutex.lock();
				if (params.is_found)
				{
					params.sbox_mutex.unlock();
					return;
				}

				sbox = params.best_sbox;
				params.sbox_mutex.unlock();

				std::swap(sbox.sbox[pos_1], sbox.sbox[pos_2]);

				sbox.cost = info.cost_function(info.cost_data.get(), sbox.sbox);
				bool is_good_nl = sbox.cost.nonlinearity >=
					info.target_properties[SBGEN_NONLINEARITY];
				while (is_good_nl)
				{
					if (info.is_log_enabled)
					{
						if (info.use_log_function && info.log_good_nl)
							info.log_good_nl_function(params, info);
						if (info.default_log_output)
						{
							std::cout << "cost=" 
							<< params.best_sbox.cost.cost 
							<< "	NL=" 
							<< params.best_sbox.cost.nonlinearity 
							<< "temperature=" 
							<< current_temperature << std::endl;
						}
					}

					if(!generator_utils::check_additional_properties(
							static_cast<properties_info_t>(info),sbox.sbox))
					{
						break;
					}

					params.sbox_mutex.lock();
					params.best_sbox = sbox;
					params.is_found = true;
					params.sbox_mutex.unlock();
					return;
				}

				params.sbox_mutex.lock();
				double cost_diff = 
					(double)(sbox.cost.cost - params.best_sbox.cost.cost);

				if (info.comparator(params.best_sbox, sbox))
				{
					params.best_sbox = sbox;
					accept_in_this_loop = true;

					if (info.is_log_enabled)
					{
						if (info.use_log_function && info.log_better_sbox)
							info.log_better_sbox_function(params, info);
						if (info.default_log_output)
						{
							std::cout << "cost=" 
							<< params.best_sbox.cost.cost 
							<< "	NL=" 
							<< params.best_sbox.cost.nonlinearity 
							<< "temperature=" 
							<< current_temperature 
							<< std::endl;
						}
					}
				} else
				{
					double u = dis(gen);
					if (u < exp(-cost_diff / current_temperature))
					{
						params.best_sbox = sbox;
						accept_in_this_loop = true;
					}
				}

				if(params.frozen_count/info.thread_count >=
					info.max_frozen_outer_loops)
				{
					params.sbox_mutex.unlock();
					return;
				}

				if (accept_in_this_loop == false)
				{
					params.frozen_count++;
					if (params.frozen_count/info.thread_count >=
						info.max_frozen_outer_loops)
					{
						break;
					}
				}
				else 
				{
					params.frozen_count = 0;
				}
				params.sbox_mutex.unlock();
			} // inner loop

			current_temperature *= info.alpha_parameter;
		} // outer loop
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
	std::optional<sbox_t> simulated_annealing (simulated_annealing_info_t<T>& info)
	{
		shared_info_t<T> thread_data;

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> distrib(0, 255);
		std::vector<std::thread> workers;


		for (int i = 0; i <= 255; i++)
		{
			thread_data.best_sbox.sbox[i] = i;
		}

		for (int i = 255; i > 0; i--)
		{
			int j = distrib(gen) % i;
			std::swap(thread_data.best_sbox.sbox[i], thread_data.best_sbox.sbox[j]);
		}

		thread_data.is_found = false;
		thread_data.frozen_count = 0;
		thread_data.iteration.store(0);
		thread_data.best_sbox.cost = info.cost_function(info.cost_data.get(),
			thread_data.best_sbox.sbox);

		for (int i = 0; i < info.thread_count; i++)
		{
			workers.push_back(
				std::thread(&simulated_annealing_thread_function<T>,
					std::ref(thread_data), std::ref(info))
			);
		}

		for (std::thread& worker : workers)
			worker.join();
		workers.clear();

		if (thread_data.is_found)
		{
			return thread_data.best_sbox.sbox;
		}

		return {};
	}

}; // namespace sbgen

#endif // _SBGEN_SIMULATED_ANNEALING_H_
