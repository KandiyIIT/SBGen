// Genetic algorithm implementation
#ifndef _SBGEN__GENETIC_H_
#define _SBGEN__GENETIC_H_

#include <thread>
#include <functional>
#include <optional>
#include <random>
#include <vector>
#include <iostream>
#include <mutex>
#include <atomic>

#include "cost_function.h"
#include "generator.h"
#include "genetic_internals.h"

namespace sbgen
{
	/**
	* @brief Genetic shared info
	* */
	template <typename T>
	struct genetic_shared_info_t
	{
		std::vector<sbox_info<T>> successors;
		std::mutex successors_mutex;

		population_t<T> population;
		std::mutex population_mutex;

		sbox_info<T> best_sbox;
		std::mutex sbox_mutex;

		std::atomic<bool> is_sbox_found;

		std::atomic<uint32_t> iteration;
	};

	/**
	* @brief Genetic method parameters
	* */
	template <typename T>
	struct genetic_info_t : public properties_info_t
	{
		using log_function_t =
			std::function<void(
				genetic_shared_info_t<T>& params,genetic_info_t<T>& info)>;
		// general parameters
		int32_t thread_count;
		int32_t mutants_per_parent;
		int32_t selection_count;
		int32_t child_per_parent;
		int32_t iterations_count;
		int32_t initial_population_count;
		int32_t crossover_count = 0;

		// flags
		bool use_crossover;
		bool is_log_enabled;
		bool delete_parents = false;
		bool use_log_function = false;
		bool default_log_output = true;
		bool log_good_nl = false;

		// functional parameters
		std::unique_ptr<cost_function_data_t> cost_data;
		cost_function_t<T> cost_function;
		genetic_comparator<T> comparator;// = comparators::less_nl<T>;
		selection_method_t<T> selection_method;
		crossover_method_t crossover_method;

		// log parameters
		uint32_t log_sbox_show_count = 10;
		uint32_t log_stride = 100;
		log_function_t log_good_nl_function;
	}; // struct genetic_info_t

	/**
	* @brief Genetic log function type
	**/
	template <typename T>
	using genetic_log_function_t =
		typename genetic_info_t<T>::log_function_t;

	/**
	* @brief Genetic mutation function
	*
	* One-thread genetic mutation function
	*
	* @param params
	*    shared data
	* @param info
	*    genetic method parameters
	* @returns 
	*   s-box with target properties in params
	*/
	template<typename T>
	void genetic_thread_function (genetic_shared_info_t<T>& params,
		genetic_info_t<T>& info) 
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> distrib(0, 255);

		while(1)
		{
			params.iteration.fetch_add(1);
			if (params.is_sbox_found.load() == true)
				break;

			params.successors_mutex.lock();
			if (params.successors.empty())
			{
				params.successors_mutex.unlock();
				break;
			}
			sbox_info<T> successor = params.successors.back();
			params.successors.pop_back();
			params.successors_mutex.unlock();

			params.population_mutex.lock();
			if (info.delete_parents == false)
				params.population.push(successor);
			params.population_mutex.unlock();

			for (int i = 0; i < info.mutants_per_parent; i++)
			{
				sbox_info<T> mutant = successor;

				int32_t pos_1 = 0;
				int32_t pos_2 = 0;
				while (pos_1 == pos_2)
				{
					pos_1 = distrib(gen);
					pos_2 = distrib(gen);
				}
				std::swap(mutant.sbox[pos_1], mutant.sbox[pos_2]);
				mutant.cost = info.cost_function(
					info.cost_data.get(), mutant.sbox);

				bool is_good_nl = mutant.cost.nonlinearity >=
					info.target_properties[SBGEN_NONLINEARITY];
				while (is_good_nl)
				{
					if (!generator_utils::check_additional_properties(
						static_cast<properties_info_t>(info), mutant.sbox))
					{
						break;
					}

					params.sbox_mutex.lock();
					params.best_sbox = mutant;
					params.is_sbox_found = true;
					if (info.is_log_enabled)
					{
						if (info.use_log_function && info.log_good_nl)
							info.log_good_nl_function(params, info);
						if (info.default_log_output)
						{
							std::cout << "found S-box with target properties"
							<< "	NL=" 
							<< params.best_sbox.cost.nonlinearity 
							<< std::endl;
							std::cout<<"SEARCH COST:"
								<<params.iteration.load()*(
									info.selection_count+info.crossover_count
								)<<std::endl;
						}
					}
					params.sbox_mutex.unlock();
					return;
				}

				params.population_mutex.lock();
				params.population.push(mutant);
				params.population_mutex.unlock();
			}
		};

		return;
	}

	/**
	* @brief Mult-thread genetic generator
	*
	* Genetic generator
	*
	* @param info
	*    genetic parameters
	* @returns 
	*   s-box with target properties 
	*/
	template<typename T>
	std::optional<sbox_t> genetic(genetic_info_t<T>& info)
	{
		genetic_shared_info_t<T> thread_data;
		population_t<T> population(info.comparator);
		population_t<T> empty_population(info.comparator);

		std::mutex population_mutex;
		std::atomic<int> population_count(info.initial_population_count);

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> distrib(0, 255);

		std::vector<std::thread> workers;

		thread_data.is_sbox_found = false; 

		for(int i=0;i< info.thread_count;i++)
		{
			workers.push_back(std::thread([&]
				{
					sbox_info<T> sbox;
					std::random_device rd;
					std::mt19937 gen(rd());
					std::uniform_int_distribution<int> distrib(0, 255);
					while (population_count.load() > 0) 
					{
						for (int i = 0; i <= 255; i++)
						{
							sbox.sbox[i] = i;
						}

						for (int i = 255; i > 0; i--) 
						{
							int j = distrib(gen) % i;
							std::swap(sbox.sbox[i], sbox.sbox[j]);
						}
						sbox.cost = info.cost_function(
							info.cost_data.get(), sbox.sbox);

						population_mutex.lock();
						population.push(std::move(sbox));
						population_mutex.unlock();
						population_count.fetch_sub(1);
					}
				}
			));
		}

		for (auto& thread : workers)
		{
			thread.join();
		}
		workers.clear();

		for (int i = 0; i < info.iterations_count; i++)
		{
			if (info.is_log_enabled &&
				i % info.log_stride == 0 && info.default_log_output)
			{
				std::cout << "Iteration " << i << std::endl;
				std::cout << "cost=" 
					<< population.top().cost.cost
					<< "	NL=" 
					<< population.top().cost.nonlinearity 
					<< std::endl;
			}
			thread_data.population = empty_population;
			info.selection_method(population, thread_data.successors,
				info.selection_count);

			// crossover can be parallelized
			if (info.use_crossover == true)
			{
				std::uniform_int_distribution<int32_t> popul_distr(0,
					thread_data.successors.size()-1);
				sbox_info<T> sb;

				for (int32_t j = 0; j < info.crossover_count; j++)
				{
					int32_t pos_1 = 0;
					int32_t pos_2 = 0;
					while (pos_1 == pos_2)
					{
						pos_1 = popul_distr(gen);
						pos_2 = popul_distr(gen);
					}

					for (int32_t k = 0; k < info.child_per_parent; k++)
					{
						sb.sbox = info.crossover_method(
							thread_data.successors[pos_1].sbox,
							thread_data.successors[pos_2].sbox);
						sb.cost = info.cost_function(
							info.cost_data.get(), sb.sbox);
						thread_data.successors.push_back(std::move(sb));
					}
				}
			}

			if (info.is_log_enabled
				&& i % info.log_stride == 0 && info.default_log_output)
			{
				std::cout<<"successors size: "
					<<thread_data.successors.size()<<std::endl;
				for (size_t z=0; z < thread_data.successors.size(); z++)
				{
					if (z > info.log_sbox_show_count)
						break;
					std::cout << "(" 
					<< thread_data.successors[z].cost.cost
					<< ", " 
					<< thread_data.successors[z].cost.nonlinearity 
					<< ")"<< std::endl;
				}
			}

			for (int j = 0;j < info.thread_count;j++)
			{
				workers.push_back(
					std::thread(
						&genetic_thread_function<T>, 
						std::ref(thread_data), 
						std::ref(info)
					)
				);
			}

			for (auto& thread : workers)
			{
				thread.join();
			}

			thread_data.successors.clear();
			workers.clear();

			if (thread_data.is_sbox_found == true)
			{
				return thread_data.best_sbox.sbox;
			}

			population = thread_data.population;
		}

		return {};
	}

}; // namespace sbgen

#endif // _SBGEN__GENETIC_H_
