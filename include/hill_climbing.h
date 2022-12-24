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

namespace sbgen
{
	/**
	* @brief Hill climbing method parameters
	**/
	template <typename T>
	struct hill_climbing_info_t : public properties_info_t
	{
		using log_function_t =
			std::function<void(shared_info_t<T>&, hill_climbing_info_t<T>&)>;
		int32_t thread_count;
		int32_t try_per_thread;
		int32_t max_frozen_count;
		bool is_log_enabled = false;
		bool use_log_function = false;
		bool default_log_output = true;
		bool log_good_nl = false;
		bool log_better_sbox = false;

		// functional parameters
		std::unique_ptr<cost_function_data_t> cost_data;
		cost_function_t<T> cost_function;
		comparator_t<T> comparator = comparators::less_nl<T>;
		log_function_t log_good_nl_function;
		log_function_t log_better_sbox_function;
	}; // struct hill_climbing_info_t

	/**
	* @brief Hill climbing log function type
	**/
	template <typename T>
	using hill_climbing_log_function_t =
		typename hill_climbing_info_t<T>::log_function_t;

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
	void hill_climbing_thread_function (shared_info_t<T>& params, 
		hill_climbing_info_t<T>& info) 
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int32_t> distrib(0, 255);

		//sbox_t new_sbox;
		//cost_info_t<T> new_cost;
		sbox_info<T> sbox;

		for (int32_t i = 0; i < info.try_per_thread; i++)
		{
			params.iteration.fetch_add(1);

			int32_t pos_1 = 0;
			int32_t pos_2 = 0;

			params.sbox_mutex.lock();
			sbox = params.best_sbox;
			params.sbox_mutex.unlock();

			while (pos_1 == pos_2)
			{
				pos_1 = distrib(gen);
				pos_2 = distrib(gen);
			}
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
						<< "	iteration=" 
						<< params.iteration.load() << std::endl;
					}
				}
			
				if (!generator_utils::check_additional_properties(
						static_cast<properties_info_t>(info),sbox.sbox))
				{
					break;
				}

				params.sbox_mutex.lock();
				params.best_sbox = sbox;
				params.is_found = true;
				if (info.is_log_enabled)
					std::cout<<"SEARCH COST:"
						<<params.iteration.load()<<std::endl;
				params.sbox_mutex.unlock();
				return;
			}

			if (!is_good_nl)
			{
				params.sbox_mutex.lock();
				params.frozen_count++;
				if (params.frozen_count > info.max_frozen_count)
				{
					if (info.is_log_enabled)
					{
						if (info.default_log_output)
						{
							std::cout << "iteration=" 
							<<params.iteration.load() 
							<< "	Search stopped: frozen_count > max_frozen_count"
							<< std::endl;
						}
					}
					params.sbox_mutex.unlock();
					return;
				}
				params.sbox_mutex.unlock();
			}

			params.sbox_mutex.lock();
			if (params.is_found)
			{
				params.sbox_mutex.unlock();
				return;
			}

			if (info.comparator(params.best_sbox, sbox))
			{
				params.best_sbox = sbox;
				params.frozen_count = 0;

				if (info.is_log_enabled)
				{
					if (info.use_log_function && info.log_better_sbox)
						info.log_better_sbox_function(params, info);
					if (info.default_log_output)
					{
						std::cout << "cost=" 
						<< params.best_sbox.cost.cost 
						<< "	NL=" << params.best_sbox.cost.nonlinearity 
						<< "	iteration=" 
						<< params.iteration.load() << std::endl;
					}
				}
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
	std::optional<sbox_t> hill_climbing(hill_climbing_info_t<T>& info)
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
				std::thread(&hill_climbing_thread_function<T>,
					std::ref(thread_data), std::ref(info))
			);
		}

		for (std::thread& worker : workers)
			worker.join();
		workers.clear();

		if (thread_data.is_found) {
			return thread_data.best_sbox.sbox;
		}

		return {};
	}

}; // namespace sbgen

#endif // _HILL_CLIMBING_H_
