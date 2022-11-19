// Internals of Genetic  algorithms
#ifndef _SBGEN_GENETIC_DETAILS_
#define _SBGEN_GENETIC_DETAILS_

#include <queue>
#include <vector>
#include <array>
#include <functional>
#include <random>

namespace sbgen {

	template <typename T>
	struct sbox_info {
		std::array<uint8_t, 256> sbox;
		cost_info_t<T> cost;
	}; // struct sbox_info
	
	template <typename T>
	using comparator_t = 
		std::function<bool(sbox_info<T>&, sbox_info<T>&)>;
        
	class comparators {
	public:
		template <typename T>
		static bool less(sbox_info<T>& a, sbox_info<T>& b) {
			return (a.cost.cost > b.cost.cost);
		}

		template <typename T>
		static bool less_nl(sbox_info<T>& a, sbox_info<T>& b) {
			return (a.cost.cost > b.cost.cost) &&
			(a.cost.nonlinearity >= b.cost.nonlinearity);
		}
	};

	template <typename T>
	struct genetic_comparator {
	public:
		comparator_t<T> comparator = comparators::less<T>;

		bool operator()(sbox_info<T>& a, sbox_info<T>& b) {
			return comparator(a,b);
		}
	}; // struct genetic_comparator

	template <typename T>
	using population_t = std::priority_queue<sbox_info<T>, 
		std::vector<sbox_info<T>>, genetic_comparator<T>>;

    template <typename T>
	using selection_method_t =
		std::function<void(population_t<T>&, std::vector<sbox_info<T>>&,int)>;

	using crossover_method_t = 
		std::function<std::array<uint8_t, 256>(
			std::array<uint8_t, 256>&, std::array<uint8_t, 256>&)>;

	class selectors {
	public:
		template <typename T>
		static void basic_selection(population_t<T>& population, 
			std::vector<sbox_info<T>>& new_population, int child_per_parent) {
			T cost_prev;

			for (int j = 0; j < child_per_parent; j++)
			{
				if (population.empty())
					return;

				sbox_info<T> s = population.top();
				population.pop();

				new_population.push_back(s);
				cost_prev = s.cost.cost;

				while(!population.empty() &&  
					population.top().cost.cost == cost_prev) {
					population.pop();
				}
			}
		}

		template <typename T>
		static void rank_selection(population_t<T>& population, 
			std::vector<sbox_info<T>>& new_population, int child_per_parent) 
		{
			int total_selected_count = 0;
			T cost_prev;
			std::random_device rd;
			std::mt19937 gen(rd());

			std::uniform_real_distribution<double> rdistrib(0.0, 1.0);
			std::vector<bool> index(population.size(),false);

			std::vector<sbox_info<T>> res;
			while (!population.empty()) {
				sbox_info<T> s = population.top();
				cost_prev = s.cost.cost;
				res.push_back(s);
				while(!population.empty() &&  
					population.top().cost.cost == cost_prev) {
					population.pop();
				}
			}
			std::uniform_int_distribution<int> distrib(0, res.size()-1);
			
			if (child_per_parent > res.size())
			{
				new_population = res;
				return;
			}

			while (total_selected_count < child_per_parent) 
			{
				int pos = distrib(gen);
				if (index[pos] == true)
					continue;

				if (rdistrib(gen) < ((2.0 * pos) / (
					child_per_parent * (child_per_parent + 1.0)))) 
				{
					index[pos] = true;
					new_population.push_back(res[res.size()-pos-1]);

					total_selected_count++;
				}
			}

			return;
		}

		template <typename T>
		static void roulette_wheel_selection(population_t<T>& population, 
			std::vector<sbox_info<T>>& new_population, int child_per_parent) 
		{
			int total_selected_count = 0;
			T Z = 0;
			T cost_prev;
			std::random_device rd;
			std::mt19937 gen(rd());

			std::uniform_real_distribution<double> rdistrib(0.0, 1.0);
			std::vector<bool> index(population.size(), false);

			std::vector<sbox_info<T>> res;
			while (!population.empty()) {
				sbox_info<T> s = population.top();
				cost_prev = s.cost.cost;
				res.push_back(s);
				while(!population.empty() &&  
					population.top().cost == cost_prev) {
					population.pop();
				}
			}
			std::uniform_int_distribution<int> distrib(0, res.size() - 1);

			if (child_per_parent > res.size())
			{
				new_population = res;
				return;
			}

			while (total_selected_count < child_per_parent)
			{
				int pos = distrib(gen);
				if (index[pos] == true)
					continue;

				if (rdistrib(gen) < ((double)res[pos].get_cost()/Z)) {
					index[pos] = true;
					new_population.push_back(res[res.size() - pos - 1]);

					total_selected_count++;
				}
			}

			return;
		}
		
	}; // class selectors

	class cossovers {
	public:
		template <typename T>
		std::array<uint8_t, 256> cycle(std::array<uint8_t, 256>& a,
										std::array<uint8_t, 256>& b)
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<int> dist(0, 255);

			std::array<bool, 256> index;
			for (int i = 0;i < 256;i++)
				index[i] = false;
			std::array<uint8_t, 256> res;

			int cycle_start = dist(gen);
			int current_pos = cycle_start;
			do  {
				res[current_pos] = a[current_pos];
				index[current_pos] = true;
				int current_val = b[current_pos];

				for (int i = 0;i < 256;i++) {
					if (a[i] == current_val) {
						current_pos = i;
						break;
					}
				}
			} while (cycle_start != current_pos);

			for (int i = 0;i < 256;i++) {
				if (index[i] == false)
					res[i] = b[i];
			}
			return std::move(res);
		}

		std::array<uint8_t, 256> pmx(std::array<uint8_t, 256>& a,
						std::array<uint8_t, 256>& b)
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<int> dist(0, 255);

			std::array<bool, 256> index;
			std::array<bool, 256> index_values;
			std::array<uint8_t, 256> res;
			for (int i = 0;i < 256;i++) {
				index[i] = false;
				index_values[i] = false;
			}

			int start_pos = 0; 
			int end_pos = 0;
			while (start_pos == end_pos) {
				start_pos = dist(gen)%256;
				end_pos = dist(gen)%256;
				if (start_pos > end_pos)
					std::swap(start_pos, end_pos);
			}

			for (int i = start_pos;i <= end_pos;i++) {
				res[i] = a[i];
				index[i] = true;
				index_values[a[i]] = true;
			}

			for (int i = 0;i < 256;i++) {
				if (index[i] == true)
					continue;
				int value = b[i];

				while (index_values[value] == true) {
					for (int j = 0;j < 256;j++) {
						if (a[j] == value) {
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
