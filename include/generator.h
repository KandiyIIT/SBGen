// Utils and macrodefinitions for s-box generators
#ifndef  _SBGEN_GENERATOR_H_
#define _SBGEN_GENERATOR_H_

#include <mutex>
#include <atomic>
#include <memory>
#include <array>

namespace sbgen {

/**
 * @brief Holds S-box properties
 *
 * Contains information about which properties to optimize at generation time 
 * and an optional random seed for the random number generator.
 **/
struct properties_info_t {
	// Array with target properties
	int32_t target_properties[SBGEN_MAX_PROPERTIES_NUMBER];
	// Flag register wich show target properties
	uint64_t properties_config;
	// Shoud we use random or precomputed seed
	bool use_random_seed;
	// Precomputed seed
	unsigned seed;

	
  /**
   * @brief Default constructor for properties_info_t
   */
	properties_info_t() : properties_config(0), use_random_seed(true) {};
	
}; // struct properties_info_t


/**
 * @brief Holds shared multithread info for generators
 *
 * Contains information about best founded s-box and it's cost function value
 **/
template <typename T>
struct shared_info_t {
	// Best S-box
	std::array<uint8_t, 256> best_sbox;
	// Cost of best S-box
	cost_info_t<T> best_cost;
	// Mutex for acces to shared_info_t fields
	std::mutex sbox_mutex;
	// Found s-box?
	bool is_found;
	// Current count of iterations without any change
	int32_t frozen_count;
	// Total count of iterations
	std::atomic<uint32_t> iteration;
}; // struct shared_info_t 

/**
 * @brief Holds utils for s-box generators
 *
 * Contains useful functions for s-box generators
 **/
class generator_utils {
public:
	
  /**
   * @brief Check properties of s-box 
   *
   * Check all target properties of s-box without nonlinearity 
   *
   * @param info
   *    target property values
   * @param sbox
   *    target s-box
   *@returns
   *   false if s-box has worse properties than target, true otherwise
   */
	static bool check_additional_properties 
	(
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
	
}; // class generator_utils 
	
}; // namespace sbgen

#endif // _SBGEN_GENERATOR_H_
