#ifndef  _SBGEN_GENERATOR_H_
#define _SBGEN_GENERATOR_H_

#include <mutex>
#include <atomic>
#include <memory>
#include <array>

namespace sbgen {

struct properties_info_t {
	int32_t target_properties[SBGEN_MAX_PROPERTIES_NUMBER];
	uint64_t properties_config;
	bool use_random_seed;
	unsigned seed;
		
	properties_info_t() : properties_config(0), use_random_seed(true) {};
};

template <typename T>
struct shared_info_t {
	std::array<uint8_t, 256> best_sbox;
	cost_info_t<T> best_cost;
	std::mutex sbox_mutex;
	bool is_found;
	int32_t frozen_count;
	std::atomic<uint32_t> iteration;
};
	
class generator_utils {
public:
	static bool check_additional_properties (
		properties_info_t info, 
		std::array<uint8_t, 256>& sbox
		) 
	{
			if(info.properties_config & SBGEN_USE_DELTA_UNIFORMITY_FLAG) 
			{
				int32_t delta_uniformity = properties::delta_uniformity(sbox);
				if (delta_uniformity >  info.target_properties[SBGEN_DELTA_UNIFORMITY]) {
					return false;
				}
			}
				
			if(info.properties_config & SBGEN_USE_ALGEBRAIC_IMMUNITY_FLAG) 
			{
				int32_t ai = properties::algebraic_immunity(sbox);
				if (ai <  info.target_properties[SBGEN_ALGEBRAIC_IMMUNITY]) {
					return false;
				} 
			}
			return true;
	}
};
	
};

#endif
