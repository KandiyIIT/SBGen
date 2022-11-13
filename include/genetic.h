// Genetic algorithm implementation
#ifndef _GENETIC_H_
#define _GENETIC_H_

#include <thread>
#include <functional>
#include <optional>
#include <random>
#include <vector>
#include <iostream>

#include "cost_function.h"
#include "generator.h"
#include "genetic_internals.h"

namespace sbgen {

/**
 * @brief Genetic method parameters
 * */
template <typename T>
struct genetic_info_t : public properties_info_t{
	// thread count
	int32_t thread_count;
	// mutants per parent
	int32_t mutants_per_parent;
	// child per parent
	int32_t child_per_parent;
	// iterations count
	int32_t iterations_count;

    // cost function and cost function data
	std::unique_ptr< cost_function_data_t> cost_data;
	std::function< cost_info_t<T>(cost_function_data_t*, std::array<uint8_t, 256>)> cost_function;
	genetic_comparator<T> comparator;
	selection_method_t<T> selection_method;
	crossover_method_t crossover_method;
}; // struct genetic_info_t 


}; // namespace sbgen

#endif // _GENETIC_H_
