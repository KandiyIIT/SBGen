// Internals of Genetic  algorithms
#ifndef _SBGEN_GENETIC_DETAILS_
#define _SBGEN_GENETIC_DETAILS_

namespace sbgen {

	template <typename T>
	struct sbox_info {
		std::array<uint8_t, 256> sbox;
		cost_info_t<T> cost;
	}; // struct sbox_info
	
	template <typename T>
	using comparator_t = 
		std::function<bool(sbox_info<T>&, sbox_info<T>&)>;

	template <typename T>
	struct genetic_comparator {
	public:
		comparator_t comparator;

		bool operator()(sbox_info<T>&, sbox_info<T>&) {
			return comparator(a,b);
		}
	}; // struct genetic_comparator

	template <typename T>
	using population_t = std::priority_queue<sbox_info<T>, 
		std::vector<sbox_info<T>>, genetic_comparator<T>>;

		template <typename T>
	using selection_method_t =
		std::function<void(population_t<T>&, std::vector<sbox_info<T>>,int)>;

	using crossover_method_t = 
		std::function<std::array<uint8_t, 256>(
			std::array<uint8_t, 256>&, std::array<uint8_t, 256>&)>;
			
	/**
	* @brief 
	*
	* Contains information about best founded s-box and it's cost function value
	**/
	template <typename T>
	struct genetic_shared_info_t {
		std::vector<sbox_info<T>> successors;
		std::mutex successors_mutex;
		population_t<T> population;
		std::mutex population_mutex;
		// Best S-box
		std::array<uint8_t, 256> best_sbox;
		// Cost of best S-box
		cost_info_t<T> best_cost;
		// Mutex for acces to shared_info_t fields
		std::mutex sbox_mutex;
		// Found s-box?
		bool is_found;
	}; // struct genetic_shared_info_t 

	class comparators {
		template <typename T>
		static bool less(sbox_info<T>& a, sbox_info<T>& b) {
			return (a.cost < b.cost);
		}

		template <typename T>
		static bool greater(sbox_info<T>& a, sbox_info<T>& b) {
			return (a.cost > b.cost);
		}

		template <typename T>
		static bool leq(sbox_info<T>& a, sbox_info<T>& b) {
			return (a.cost <= b.cost);
		}

		template <typename T>
		static bool geq(sbox_info<T>& a, sbox_info<T>& b) {
			return (a.cost >= b.cost);
		}

		template <typename T>
		static bool less_nl(sbox_info<T>& a, sbox_info<T>& b) {
			return (a.cost < b.cost) && (a.nonlinearity <= b.nonlinearity);
		}
	};

	class selectors {
		template <typename T>
		static void basic_selection(population_t& population, 
			std::vector<sbox_info<T>>& new_population, int child_per_parent) {
			T cost_prev;

			for (int j = 0; j < child_per_parent; j++)
			{
				if (population.empty())
					return;

				sbox_info<T> s = population.top();
				population.pop();

				new_population.push_back(s);
				cost_prev = s.cost;

				while(!population.empty() &&  
					population.top().cost == cost_prev) {
					population.pop();
				}
			}
		}

		template <typename T>
		static void rank_selection(population_t& population, 
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
				cost_prev = s.cost;
				res.push_back(s);
				while(!population.empty() &&  
					population.top().cost == cost_prev) {
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
		static void roulette_wheel_selection(population_t& population, 
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
				cost_prev = s.cost;
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
				res[current_pos] = parent1[current_pos];
				index[current_pos] = true;
				int current_val = parent2[current_pos];

				for (int i = 0;i < 256;i++) {
					if (parent1[i] == current_val) {
						current_pos = i;
						break;
					}
				}
			} while (cycle_start != current_pos);

			for (int i = 0;i < 256;i++) {
				if (index[i] == false)
					res[i] = parent2[i];
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
			std::array<uint8_t, 256> res(false);
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
				res[i] = parent1[i];
				index[i] = true;
				index_values[parent1[i]] = true;
			}

			for (int i = 0;i < 256;i++) {
				if (index[i] == true)
					continue;
				int value = parent2[i];

				while (index_values[value] == true) {
					for (int j = 0;j < 256;j++) {
						if (parent1[j] == value) {
							value = parent2[j];
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
