// Internals of Genetic algorithms
#ifndef _SBGEN_GENETIC_DETAILS_
#define _SBGEN_GENETIC_DETAILS_

#include <queue>
#include <vector>
#include <array>
#include <functional>
#include <random>

#include "cost_function.h"
#include "sbox_properties.h"

namespace sbgen 
{
	/**
	 * @brief info about s-box in population
	 * */
	template <typename T>
	struct sbox_info {
		sbox_t sbox;
		cost_info_t<T> cost;
	}; // struct sbox_info

	/**
	 * @brief comparator function definition for population
	 * */
	template <typename T>
	using comparator_t = 
		std::function<bool(sbox_info<T>&, sbox_info<T>&)>;

	/**
	 * @brief comparators
	 * */
	class comparators {
	public:
		/**
		* @brief less - ordering from smallest to largest
		* @param a - first s-box
		* @param b - second s-box
		* @return true - need exchange s-boxes, other - false
		* */
		template <typename T>
		static bool less(sbox_info<T>& a, sbox_info<T>& b) {
			return (a.cost.cost > b.cost.cost);
		}

		/**
		* @brief less_nl - ordering from smallest to largest
		* taking into account the nonlinearity
		* @param a - first s-box
		* @param b - second s-box
		* @return true - need exchange s-boxes, other - false
		* */
		template <typename T>
		static bool less_nl(sbox_info<T>& a, sbox_info<T>& b) {
			if (a.cost.nonlinearity > b.cost.nonlinearity)
				return false;
			if (a.cost.nonlinearity < b.cost.nonlinearity)
				return true;
			return (a.cost.cost > b.cost.cost);
		}
	}; // class comparators

	/**
	 * @brief comparator wrapper for std::priority_queue
	 * */
	template <typename T>
	struct genetic_comparator {
	public:
		comparator_t<T> comparator = comparators::less<T>;

		bool operator()(sbox_info<T>& a, sbox_info<T>& b) {
			return comparator(a,b);
		}
	}; // struct genetic_comparator

	/**
	 * @brief population type pseudonim
	 * */
	template <typename T>
	using population_t = std::priority_queue<sbox_info<T>,
		std::vector<sbox_info<T>>, genetic_comparator<T>>;

	/**
	 * @brief selection method pseudonim
	 * */
	template <typename T>
	using selection_method_t =
		std::function<void(population_t<T>&, std::vector<sbox_info<T>>&,int)>;

	/**
	 * @brief crossover method pseudonim
	 * */
	using crossover_method_t = std::function<sbox_t(sbox_t&, sbox_t&)>;

	/**
	 * @brief selectors
	 * */
	class selectors {
	public:
		/**
		 * @brief basic selection (selects only best s-boxes)
		 * @param population - target population
		 * @param successors - selected s-boxes
		 * @param count      - target successors count
		 * @return returned value via ref in successors
		 */
		template <typename T>
		static void basic_selection(population_t<T>& population,
			std::vector<sbox_info<T>>& successors, const size_t count)
		{
			T cost_prev{};

			if (population.size() <= count)
			{
				while (!population.empty())
				{
					successors.push_back(std::move(
						const_cast<sbox_info<T>&>(population.top())));
					population.pop();
				}
				return;
			}

			for (size_t j = 0; j < count; j++)
			{
				if (population.empty())
					return;

				successors.push_back(std::move(
					const_cast<sbox_info<T>&>(population.top())));
				population.pop();
				cost_prev = successors[j].cost.cost;

				while(!population.empty() &&
					population.top().cost.cost == cost_prev)
				{
					population.pop();
				}
			}

			return;
		}

		/**
		 * @brief rank selection (selects s-boxes via rank method)
		 * @param population - target population
		 * @param successors - selected s-boxes
		 * @param count      - target successors count
		 * @return returned value via ref in successors
		 * @warning legacy code. Only for compatibility
		 */
		template <typename T>
		static void rank_selection(population_t<T>& population,
			std::vector<sbox_info<T>>& successors, const size_t count)
		{
			size_t                                 total_selected_count = 0;
			T                                      cost_prev{};
			std::random_device                     rd;
			std::mt19937                           gen(rd());
			std::uniform_real_distribution<double> rdistrib(0.0, 1.0);
			std::vector<bool>                      index(population.size(), false);
			std::vector<sbox_info<T>>              res;

			while (!population.empty())
			{
				res.push_back(std::move(
					const_cast<sbox_info<T>&>(population.top())));
				population.pop();
				cost_prev = res.back().cost.cost;

				while(!population.empty() &&  
					population.top().cost.cost == cost_prev)
				{
					population.pop();
				}
			}

			if (res.size() <= count)
			{
				successors = std::move(res);
				return;
			}

			std::uniform_int_distribution<size_t> distrib(0, res.size()-1);

			while (total_selected_count < count) 
			{
				size_t pos = distrib(gen);
				if (index[pos] == true)
					continue;

				if (rdistrib(gen) < ((2.0 * (pos)) / (count * (count + 1.0))))
				{
					index[pos] = true;
					successors.push_back(std::move(res[pos]));
					total_selected_count++;
				}
			}

			return;
		}

		/**
		 * @brief roulette wheel selection (selects s-boxes
		 * via roulette wheel method)
		 * @param population - target population
		 * @param successors - selected s-boxes
		 * @param count      - target successors count
		 * @return returned value via ref in successors
		 * @warning legacy code. Only for compatibility
		 */
		template <typename T>
		static void roulette_wheel_selection(population_t<T>& population,
			std::vector<sbox_info<T>>& successors, const size_t count)
		{
			size_t                                 total_selected_count = 0;
			T                                      cost_sum{};
			T                                      cost_prev{};
			std::random_device                     rd;
			std::mt19937                           gen(rd());
			std::uniform_real_distribution<double> rdistrib(0.0, 1.0);
			std::vector<bool>                      index(population.size(), false);
			std::vector<sbox_info<T>>              res;

			while (!population.empty())
			{
				res.push_back(std::move(
					const_cast<sbox_info<T>&>(population.top())));
				population.pop();
				cost_prev = res.back().cost.cost;
				cost_sum += cost_prev;

				while(!population.empty() &&  
					population.top().cost.cost == cost_prev) 
				{
					population.pop();
				}
			}

			if (res.size() <= count)
			{
				successors = std::move(res);
				return;
			}

			std::uniform_int_distribution<size_t> distrib(0, res.size()-1);

			while (total_selected_count < count)
			{
				size_t pos = distrib(gen);
				if (index[pos] == true)
					continue;

				if (rdistrib(gen) < ((double)res[pos].cost.cost/cost_sum))
				{
					index[pos] = true;
					successors.push_back(std::move(res[pos]));
					total_selected_count++;
				}
			}

			return;
		}

		/**
		 * @brief modified rank selection (selects s-boxes via rank method)
		 * @param population - target population
		 * @param successors - selected s-boxes
		 * @param count      - target successors count
		 * @return returned value via ref in successors
		 */
		template <typename T>
		static void rank_sequential_selection(population_t<T>& population,
			std::vector<sbox_info<T>>& successors, const size_t count)
		{
			size_t                                 total_selected_count = 0;
			size_t                                 current_position = 0;
			T                                      cost_prev{};
			std::random_device                     rd;
			std::mt19937                           gen(rd());
			std::uniform_real_distribution<double> rdistrib(0.0, 1.0);
			std::vector<bool>                      index(population.size(),false);
			std::vector<sbox_info<T>>              res;

			while (!population.empty())
			{
				res.push_back(std::move(
					const_cast<sbox_info<T>&>(population.top())));
				population.pop();
				cost_prev = res.back().cost.cost;

				while(!population.empty() &&
					population.top().cost.cost == cost_prev)
				{
					population.pop();
				}
			}

			if (res.size() <= count)
			{
				successors = std::move(res);
				return;
			}

			std::uniform_int_distribution<size_t> distrib(0, res.size()-1);

			while (total_selected_count < count) 
			{
				if (index[current_position] == true)
					continue;

				if (rdistrib(gen) < (1 - (2.0 * (current_position)) / (
					count * (count + 1.0))))
				{
					index[current_position] = true;
					successors.push_back(std::move(res[current_position]));
					total_selected_count++;
				}
				current_position = current_position + 1 % res.size();
			}

			return;
		}

		/**
		 * @brief modified roulette wheel selection
		 * (selects s-boxes via roulette wheel method)
		 * @param population - target population
		 * @param successors - selected s-boxes
		 * @param count      - target successors count
		 * @return returned value via ref in successors
		 */
		template <typename T>
		static void roulette_wheel_sequential_selection(
			population_t<T>& population, std::vector<sbox_info<T>>& successors,
			const size_t count) 
		{
			size_t                                 total_selected_count = 0;
			size_t                                 current_position = 0;
			T                                      cost_sum{};
			T                                      cost_prev{};
			std::random_device                     rd;
			std::mt19937                           gen(rd());
			std::uniform_real_distribution<double> rdistrib(0.0, 1.0);
			std::vector<bool>                      index(population.size(), false);
			std::vector<sbox_info<T>>              res;

			while (!population.empty())
			{
				res.push_back(std::move(
					const_cast<sbox_info<T>&>(population.top())));
				population.pop();
				cost_prev = res.back().cost.cost;
				cost_sum += cost_prev;

				while(!population.empty() &&
					population.top().cost.cost == cost_prev)
				{
					population.pop();
				}
			}

			if (res.size() <= count)
			{
				successors = std::move(res);
				return;
			}

			std::uniform_int_distribution<size_t> distrib(0, res.size() - 1);

			while (total_selected_count < count)
			{
				if (index[current_position] == true)
					continue;

				if (rdistrib(gen) < ((double)(
					1 - res[current_position].cost.cost/cost_sum)))
				{
					index[current_position] = true;
					successors.push_back(std::move(res[current_position]));
					total_selected_count++;
				}
				current_position = current_position + 1 % res.size();
			}

			return;
		}
	}; // class selectors

	/**
	 * @brief cossovers
	 * */
	class cossovers {
	public:
		/**
		 * @brief cycle crossover
		 * @param a - first parent
		 * @param b - second parent
		 * @return "child" s-box
		 * */
		static sbox_t cycle(sbox_t& a, sbox_t& b)
		{
			std::random_device                 rd;
			std::mt19937                       gen(rd());
			std::uniform_int_distribution<int> dist(0, 255);
			std::vector<bool>                  index(256, false);
			sbox_t                             res;
			int32_t                            cycle_start = dist(gen);
			int32_t                            current_pos = cycle_start;

			do
			{
				res[current_pos] = a[current_pos];
				index[current_pos] = true;
				int current_val = b[current_pos];

				for (int i = 0; i < 256; i++)
				{
					if (a[i] == current_val)
					{
						current_pos = i;
						break;
					}
				}
			} while (cycle_start != current_pos);

			for (int i = 0; i < 256; i++)
			{
				if (index[i] == false)
					res[i] = b[i];
			}

			return std::move(res);
		}

		/**
		 * @brief pmx crossover
		 * @param a - first parent
		 * @param b - second parent
		 * @return "child" s-box
		 * */
		static sbox_t pmx(sbox_t& a, sbox_t& b)
		{
			std::random_device                 rd;
			std::mt19937                       gen(rd());
			std::uniform_int_distribution<int> dist(0, 255);
			std::vector<bool>                  index(256, false);
			std::vector<bool>                  index_values(256, false);
			sbox_t                             res;
			int32_t                            start_pos = 0;
			int32_t                            end_pos = 0;

			while (start_pos == end_pos)
			{
				start_pos = dist(gen);
				end_pos = dist(gen);
				if (start_pos > end_pos)
					std::swap(start_pos, end_pos);
			}

			for (int i = start_pos; i <= end_pos; i++)
			{
				res[i] = a[i];
				index[i] = true;
				index_values[a[i]] = true;
			}

			for (int i = 0;i < 256;i++)
			{
				if (index[i] == true)
					continue;

				uint8_t value = b[i];

				while (index_values[value] == true)
				{
					for (int j = 0;j < 256;j++)
					{
						if (a[j] == value) 
						{
							value = b[j];
							break;
						}
					}
				}

				res[i] = value;
				index[i] = true;
				index_values[value] = true;
			}

			return std::move(res);
		}
	}; // class cossovers

}; // namespace sbgen

#endif // _SBGEN_GENETIC_DETAILS_
