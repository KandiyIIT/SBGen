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
	* @brief Genetic method parameters
	* */
	template <typename T>
	struct genetic_info_t : public properties_info_t
	{
		int32_t thread_count;
		int32_t mutants_per_parent;
		int32_t selection_count;
		int32_t child_per_parent;
		int32_t iterations_count;
		int32_t initial_population_count;
		int32_t crossover_count;
		bool use_crossover;
		bool is_log_enabled;
		bool delete_parents;

		// cost function and cost function data
		std::unique_ptr<cost_function_data_t> cost_data;
		cost_function_t cost_function;
		genetic_comparator<T> comparator;
		selection_method_t<T> selection_method;
		crossover_method_t crossover_method;

		// log settings
		uint32_t log_sbox_show_count;
		uint32_t log_print_best_sbox;
		uint32_t log_stride;
	}; // struct genetic_info_t 

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
	};

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
	void genetic_thread_function (
		genetic_shared_info_t<T>& params, 
		genetic_info_t<T>& info
	) 
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> distrib(0, 255);
		//std::array<uint8_t, 256> new_sbox;

		while(1)
		{
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
				auto new_cost = info.cost_function(
					info.cost_data.get(), mutant.sbox);
				mutant.cost = new_cost;

				bool is_good_nl = new_cost.nonlinearity >=
					info.target_properties[SBGEN_NONLINEARITY];
				while (is_good_nl)
				{
					if (info.is_log_enabled)
					{
						std::cout << "new best cost=" 
						<< params.best_sbox.cost.cost
						<< "	NL=" 
						<< params.best_sbox.cost.nonlinearity 
						<< std::endl;
					}

					if(!generator_utils::check_additional_properties(
						static_cast<properties_info_t>(info),mutant.sbox))
					{
						break;
					}
					params.sbox_mutex.lock();
					params.best_sbox = mutant;
					params.best_sbox.cost = new_cost;
					params.is_sbox_found = true;
					if (info.is_log_enabled)
					{
						std::cout << "found"
						<< "	NL=" 
						<< params.best_sbox.cost.nonlinearity 
						<< std::endl;
					}
					params.sbox_mutex.unlock();
					return;
				}
				params.population_mutex.lock();
				params.population.push(mutant);
				params.population_mutex.unlock();
			} // for (i)
		}; // while(1)

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
	std::optional<std::array<uint8_t, 256>> genetic(genetic_info_t<T>& info)
	{
		genetic_shared_info_t<T> thread_data;
		std::random_device rd;

		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> distrib(0, 255);
		std::vector<std::thread> workers;
		population_t<T> population(info.comparator);
		population_t<T> empty_population(info.comparator);
		std::mutex population_mutex;
		std::atomic<int> population_count(info.initial_population_count);
		thread_data.is_sbox_found = false; 

		for(int i=0;i< info.thread_count;i++) {
			workers.push_back(std::thread([&]
				{
					sbox_info<T> sbox;
					std::random_device rd;
					std::mt19937 gen(rd());
					std::uniform_int_distribution<int> distrib(0, 255);
					while (population_count.load() > 0) 
					{
						population_mutex.lock();
						for (int i = 0; i <= 255; i++) {
							sbox.sbox[i] = i;
						}

						for (int i = 255; i > 0; i--) 
						{
							int j = distrib(gen) % i;
							std::swap(sbox.sbox[i], sbox.sbox[j]);
						}
						if (properties::is_bijective(sbox.sbox) == false)
						{
							printf("ERROR 3\n");
							for(;;);
						}

						sbox.cost = info.cost_function(
							info.cost_data.get(), sbox.sbox);
						population.push(sbox);
						population_mutex.unlock();
						population_count.fetch_sub(1);
					}
				}));
		}

		for (auto& thread : workers) {
			thread.join();
		}

		workers.clear();

		for (int i = 0; i < info.iterations_count; i++)
		{
			thread_data.best_sbox = population.top();
			if (info.is_log_enabled)
			{
				if (i % 100 == 0)
				{
					std::cout << "Iteration " << i << std::endl;
					if (info.is_log_enabled)
					{
						std::cout << "cost=" 
						<< thread_data.best_sbox.cost.cost
						<< "	NL=" 
						<< thread_data.best_sbox.cost.nonlinearity 
						<< std::endl;
					}
				}
			}
			thread_data.population = empty_population;
			info.selection_method(population, thread_data.successors,
				info.selection_count);

			// TODO paralelize crossover
			if (info.use_crossover == true)
			{
				std::uniform_int_distribution<int> popul_distr(0,
					thread_data.successors.size()-1);
				sbox_info<T> sb;

				for (int j = 0; j < info.crossover_count; j++)
				{
					int pos_1 = 0;
					int pos_2 = 0;

					while (pos_1 == pos_2)
					{
						pos_1 = popul_distr(gen);
						pos_2 = popul_distr(gen);
					}

					for (int k = 0; k < info.child_per_parent; k++)
					{
						if (properties::is_bijective(
							thread_data.successors[pos_1].sbox) == false ||
						properties::is_bijective(
							thread_data.successors[pos_2].sbox) == false)
						{
							printf("ERROR 1!\n");
							for(;;);
						}
						sb.sbox = info.crossover_method(
							thread_data.successors[pos_1].sbox,
							thread_data.successors[pos_2].sbox);
						sb.cost = info.cost_function(
							info.cost_data.get(), sb.sbox);
						thread_data.successors.push_back(sb);
					}
				}
			}//crossover

			if (info.is_log_enabled)
			{
				if (i % 100 == 0)
				{
					std::cout<<"successors size: "
						<<thread_data.successors.size()<<std::endl;
					for (size_t z=0; z < thread_data.successors.size(); z++)
					{
						if (z > 10)
							break;
					std::cout << "(" 
					<< thread_data.successors[z].cost.cost
					<< ", " 
					<< thread_data.successors[z].cost.nonlinearity 
					<< ")"<< std::endl;
					}
				}
			}

			for (int j = 0;j < info.thread_count;j++) {
				workers.push_back(
					std::thread(
						&genetic_thread_function<T>, 
						std::ref(thread_data), 
						std::ref(info)
					)
				);
			}

			for (auto& thread : workers) {
				thread.join();
			}

			thread_data.successors.clear();
			workers.clear();

			if (thread_data.is_sbox_found == true) {
				return thread_data.best_sbox.sbox;
			}

			population = thread_data.population;
		}

		return {};
	}

}; // namespace sbgen

#endif // _SBGEN__GENETIC_H_
