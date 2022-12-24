// Cost functions implementation
#ifndef _COST_FUNCTION_H_
#define _COST_FUNCTION_H_

#include <algorithm>
#include <array>
#include <cstdint>
#include <cassert>

#include "sbox_properties.h"

namespace sbgen
{

	/**
	* @brief Holds S-box cost and nonlinearity
	**/
	template <typename T>
	struct cost_info_t {
		T cost;
		int32_t nonlinearity;
	};

	/**
	 * @brief info about s-box
	 * */
	template <typename T>
	struct sbox_info {
		sbox_t sbox;
		cost_info_t<T> cost;
	}; // struct sbox_info

	/**
	 * @brief comparator function definition
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
	* @brief Base class for cost functions parameters
	**/
	class cost_function_data_t {
	public:
		/**
		* @brief Cost function name
		**/
		virtual std::string name() {
				return "none";
		}
	};

	/**
	 * @brief cost function type 
	 **/
	template <typename T>
	using cost_function_t = std::function< cost_info_t<T>(
		cost_function_data_t*, sbox_t)>;

	/**
	* @brief WHS cost function parameters
	**/
	class whs_function_data_t : public cost_function_data_t
	{
	public:
		int32_t r; // r-parameter 
		int32_t x; // x-parameter

		/**
		* @brief Cost function name
		**/
		std::string name() override{
				return "whs";
		}

		/**
		* @brief Constructor for whs_function_data_t
		*/
		whs_function_data_t(int _r, int _x) : r(_r), x(_x) {};

		/**
		* @brief Destructor for whs_function_data_t
		*/
		virtual ~whs_function_data_t() {};

	}; // class whs_function_data_t


	/**
	* @brief WHS cost function
	*
	* Sum( abs(abs(WHT(i,j)) - x)^r) ) for i = 0..255 , j= 0..7    
	*
	* @param _data
	*    whs cost function data
	* @param sbox
	*    target s-box
	*@returns
	*   cost of target s-box
	*/
	// specialization for case r<0
	template <typename T>
	cost_info_t<T> _whs_minus (
		cost_function_data_t* _data, 
		sbox_t sbox
	) 
	{
		uint8_t truth_table[256];
		int32_t max_spectre = 0;
		int spectre[256];
		cost_info_t<T> cost;

		whs_function_data_t* data = static_cast<whs_function_data_t*>(_data);
		if (data->r < 0)
			data->r = -data->r;

		cost.cost = 0;
		for (int b = 1; b < 256; b++)
		{
			for (int i = 0; i < 256; i++)
			{
				truth_table[i] =
					sbgen::transform_utils::one_bits[sbox[i] & b] & 0x01;
			}

			sbgen::transform_utils::fwht_transform(truth_table, spectre);

			for (int i = 0; i < 256; i++)
			{
				int spectre_temp = spectre[i];

				if (spectre_temp < 0)
					spectre_temp = -spectre_temp;
				T val = ((spectre_temp - data->x) >= 0) ?
					(spectre_temp - data->x) : -(spectre_temp - data->x);

				assert(((void)"error: absolute spectre value less 0", val >= 0));

				T part = val;
				T one = 1.0;

				for (int k = 1; k < data->r; k++)
					part *= val;
				if (part != 0)
					cost.cost += one/part;

				if (spectre_temp > max_spectre)
					max_spectre = static_cast<int32_t>(spectre_temp);
			}
		}

		cost.nonlinearity = 128 - max_spectre / 2;

		return cost;
	}

	template <typename T>
	cost_info_t<T> whs (
		cost_function_data_t* _data, 
		sbox_t sbox
	) 
	{
		uint8_t truth_table[256];
		int32_t max_spectre = 0;
		int spectre[256];
		cost_info_t<T> cost;

		whs_function_data_t* data = static_cast<whs_function_data_t*>(_data);
		if (data->r < 0)
			return _whs_minus<T>(_data, sbox);

		cost.cost = 0;
		for (int b = 1; b < 256; b++)
		{
			for (int i = 0; i < 256; i++)
			{
				truth_table[i] =
					sbgen::transform_utils::one_bits[sbox[i] & b] & 0x01;
			}

			sbgen::transform_utils::fwht_transform(truth_table, spectre);

			for (int i = 0; i < 256; i++)
			{
				int spectre_temp = spectre[i];

				if (spectre_temp < 0)
					spectre_temp = -spectre_temp;
				T val = ((spectre_temp - data->x) >= 0) ?
					(spectre_temp - data->x) : -(spectre_temp - data->x);

				assert(((void)"error: absolute spectre value less 0", val >= 0));

				T part = val;

				for (int k = 1; k < data->r; k++)
					part *= val;
				cost.cost += part;

				if (spectre_temp > max_spectre)
					max_spectre = static_cast<int32_t>(spectre_temp);
			}
		}

		cost.nonlinearity = 128 - max_spectre / 2;

		return cost;
	}

	/**
	* @brief WCF cost function parameters
	**/
	class wcf_function_data_t : public cost_function_data_t {
	public:
		/**
		* @brief Cost function name
		**/
		std::string name() override{
				return "wcf";
		}

		/**
		* @brief Destructor for wcf_function_data_t
		*/
		virtual ~wcf_function_data_t() {};
	};

	/**
	* @brief WCF cost function
	*
	* Sum( abs(WHT(i,j)*(WHT(i,j)-4)*...*(WHT(i,j)-32)) ) for i = 0..255 , j= 0..7    
	*
	* @param _data
	*    whs cost function data
	* @param sbox
	*    target s-box
	*@returns
	*   cost of target s-box
	*/
	template <typename T>
	cost_info_t<T> wcf (
		cost_function_data_t* _data, 
		sbox_t sbox
	) 
	{
		uint8_t truth_table[256];
		int32_t max_spectre = 0;
		int spectre[256];
		cost_info_t<T> cost;

		cost.cost = 0;
		for (int b = 1; b < 256; b++)
		{
			for (int i = 0; i < 256; i++)
			{
				truth_table[i] =
					sbgen::transform_utils::one_bits[sbox[i] & b] & 0x01;
			}

			sbgen::transform_utils::fwht_transform(truth_table, spectre);

			for (int i = 0; i < 256; i++)
			{
				int spectre_temp = spectre[i];

				if (spectre_temp < 0)
					spectre_temp = -spectre_temp;
                
				if(spectre_temp<=32)
					continue;

				T part = 1;

				for (int k = 32; k >= 0; k-=4)
					part *= (spectre_temp-k);
				cost.cost += part;

				if (spectre_temp > max_spectre)
					max_spectre = static_cast<int32_t>(spectre_temp);
			}
		}

		cost.nonlinearity = 128 - max_spectre / 2;

		return cost;
	}

	/**
	* @brief PCF cost function parameters
	**/
	class pcf_function_data_t : public cost_function_data_t {
	public:
		int32_t level;

		/**
		* @brief Cost function name
		**/
		std::string name() override{
				return "pcf";
		}
        
        
		/**
		* @brief Constructor for pcf_function_data_t
		*/
		pcf_function_data_t(int32_t _level) : level(_level) {};

		/**
		* @brief Destructor for pcf_function_data_t
		*/
		virtual ~pcf_function_data_t() {};
	};

	/**
	* @brief PCF cost function
	*
	* Sum( 2^-i * H(i-l) ) where i=0...N , H - historgam of absolute WHT , l - largest nonzero in WHT
	*
	* @param _data
	*    whs cost function data
	* @param sbox
	*    target s-box
	*@returns
	*   cost of target s-box
	*/
	template <typename T>
	cost_info_t<T> pcf (
		cost_function_data_t* _data, 
		sbox_t sbox
	) 
	{
		uint8_t					truth_table[256];
		int32_t					max_spectre = 0;
		int						spectre[256];
		int32_t					Histogram[256+1]={0,};
		uint32_t 				max_index = 256;
		cost_info_t<T> 	cost;
        
		pcf_function_data_t* data = static_cast<pcf_function_data_t*>(_data);

		cost.cost = 0;
		for (int b = 1; b < 256; b++)
		{
			for (int i = 0; i < 256; i++)
			{
				truth_table[i] =
					sbgen::transform_utils::one_bits[sbox[i] & b] & 0x01;
			}


			sbgen::transform_utils::fwht_transform(truth_table, spectre);

			for (int i = 0; i < 256; i++)
			{
				int spectre_temp = spectre[i];

				if (spectre_temp < 0)
					spectre_temp = -spectre_temp;
                
				Histogram[spectre_temp]++;

				if (spectre_temp > max_spectre)
					max_spectre = static_cast<int32_t>(spectre_temp);
			}
		}

		while(Histogram[max_index]==0)
			max_index-=4;

		for(int i=0; i < data->level; i++)
			cost.cost+=((double)Histogram[max_index - i])/(1<<i);

		cost.nonlinearity = 128 - max_spectre / 2;

		return cost;
	}

	/**
	* @brief CF1 cost function parameters
	**/
	class cf1_function_data_t : public cost_function_data_t {
	public:
		int32_t r; // r-parameter 
		int32_t x; // x-parameter
		int32_t y; // x-parameter

		/**
		* @brief Cost function name
		**/
		std::string name() override{
				return "cf1";
		}

		/**
		* @brief Constructor for cf1_function_data_t
		*/
		cf1_function_data_t(int _r, int _x, int _y) : r(_r), x(_x), y(_y) {};

		/**
		* @brief Destructor for cf1_function_data_t
		*/
		virtual ~cf1_function_data_t() {};
	}; // class cf1_function_data_t


	/**
	* @brief CF1 cost function
	*
	* Sum( abs(abs(WHT(i,j)) - y)^r)/4 , [[WHT(i,j) < X]] )
	* for i = 0..255 , j= 0..7   
	*
	* @param _data
	*    cf1 cost function data
	* @param sbox
	*    target s-box
	*@returns
	*   cost of target s-box
	*/
	// specialization for case r < 0
	template <typename T>
	cost_info_t<T> _cf1_minus (
		cost_function_data_t* _data, 
		sbox_t sbox
	) 
	{
		uint8_t truth_table[256];
		int32_t max_spectre = 0;
		int spectre[256];
		cost_info_t<T> cost;

		cf1_function_data_t* data = static_cast<cf1_function_data_t*>(_data);
		if (data->r < 0)
			data->r = -data->r;

		cost.cost = 0;
		for (int b = 1; b < 256; b++)
		{
			for (int i = 0; i < 256; i++)
			{
				truth_table[i] =
					sbgen::transform_utils::one_bits[sbox[i] & b] & 0x01;
			}


			sbgen::transform_utils::fwht_transform(truth_table, spectre);

			for (int i = 0; i < 256; i++)
			{
				int spectre_temp = spectre[i];

				if (spectre_temp < 0)
					spectre_temp = -spectre_temp;

				if (spectre_temp <= data->x)
					continue;

				T val = ((spectre_temp - data->y) >= 0) ? 
					(spectre_temp - data->y) : -(spectre_temp - data->y);

				assert(((void)"error: absolute spectre value less 0", val >= 0));

				T part = val;
				T one = 1.0;

				for (int k = 1; k < data->r; k++)
					part *= val;
				if (part != 0)
					cost.cost += (one/part)*4;

				if (spectre_temp > max_spectre)
					max_spectre = static_cast<int32_t>(spectre_temp);
			}
		}

		cost.nonlinearity = 128 - max_spectre / 2;

		return cost;
	}

	template <typename T>
	cost_info_t<T> cf1 (
		cost_function_data_t* _data, 
		sbox_t sbox
	) 
	{
		uint8_t truth_table[256];
		int32_t max_spectre = 0;
		int spectre[256];
		cost_info_t<T> cost;

		cf1_function_data_t* data = static_cast<cf1_function_data_t*>(_data);
		if (data->r < 0)
			return _cf1_minus<T>(_data, sbox);

		cost.cost = 0;
		for (int b = 1; b < 256; b++)
		{
			for (int i = 0; i < 256; i++)
			{
				truth_table[i] =
					sbgen::transform_utils::one_bits[sbox[i] & b] & 0x01;
			}


			sbgen::transform_utils::fwht_transform(truth_table, spectre);

			for (int i = 0; i < 256; i++)
			{
				int spectre_temp = spectre[i];

				if (spectre_temp < 0)
					spectre_temp = -spectre_temp;

				if (spectre_temp <= data->x)
					continue;

				T val = ((spectre_temp - data->y) >= 0) ? 
					(spectre_temp - data->y) : -(spectre_temp - data->y);

				assert(((void)"error: absolute spectre value less 0", val >= 0));

				T part = val;

				for (int k = 1; k < data->r; k++)
					part *= val;
				cost.cost += part/4;

				if (spectre_temp > max_spectre)
					max_spectre = static_cast<int32_t>(spectre_temp);
			}
		}

		cost.nonlinearity = 128 - max_spectre / 2;

		return cost;
	}

	/**
	* @brief CF2 cost function parameters
	**/
	class cf2_function_data_t : public cost_function_data_t {
	public:
		int32_t r; // r-parameter 
		int32_t x; // x-parameter
		int32_t y; // x-parameter

		/**
		* @brief Cost function name
		**/
		std::string name() override{
				return "cf2";
		}

		/**
		* @brief Constructor for cf1_function_data_t
		*/
		cf2_function_data_t(int _r, int _x, int _y) : r(_r), x(_x), y(_y) {};
	
		/**
		* @brief Destructor for cf1_function_data_t
		*/
		virtual ~cf2_function_data_t() {};
	}; // class cf_function_data_t


	/**
	* @brief CF2 cost function
	*
	* Sum( 2^(abs(abs(WHT(i,j)) - y)*R)/4) , [[WHT(i,j) < X]] )
	* for i = 0..255 , j= 0..7   
	*
	* @param _data
	*    cf2 cost function data
	* @param sbox
	*    target s-box
	*@returns
	*   cost of target s-box
	*/
	template <typename T>
	cost_info_t<T> cf2 (
		cost_function_data_t* _data, 
		sbox_t sbox
	) 
	{
		uint8_t truth_table[256];
		int32_t max_spectre = 0;
		int spectre[256];
		cost_info_t<T> cost;

		cf2_function_data_t* data = static_cast<cf2_function_data_t*>(_data);

		cost.cost = 0;
		for (int b = 1; b < 256; b++)
		{
			for (int i = 0; i < 256; i++)
			{
				truth_table[i] =
					sbgen::transform_utils::one_bits[sbox[i] & b] & 0x01;
			}

			sbgen::transform_utils::fwht_transform(truth_table, spectre);

			for (int i = 0; i < 256; i++)
			{
				int spectre_temp = spectre[i];

				if (spectre_temp < 0)
					spectre_temp = -spectre_temp;

				if (spectre_temp <= data->x)
					continue;

				int32_t val = ((spectre_temp - data->y) >= 0) ? 
					(spectre_temp - data->y) : -(spectre_temp - data->y);

				uint64_t part = (val>>2)*data->r;

				assert(((void)"error: type overflow", part < 64));
				assert(((void)"error: negative value", part >= 0));

				cost.cost += static_cast<T>(1<<part);

				if (spectre_temp > max_spectre)
					max_spectre = static_cast<int32_t>(spectre_temp);
			}
		}

		cost.nonlinearity = 128 - max_spectre / 2;

		return cost;
	}

	template <>
	cost_info_t<double> cf2 (
		cost_function_data_t* _data, 
		sbox_t sbox
	) 
	{
		uint8_t truth_table[256];
		int32_t max_spectre = 0;
		int spectre[256];
		cost_info_t<double> cost;

		cf2_function_data_t* data = static_cast<cf2_function_data_t*>(_data);

		cost.cost = 0;
		for (int b = 1; b < 256; b++)
		{
			for (int i = 0; i < 256; i++)
			{
				truth_table[i] =
					sbgen::transform_utils::one_bits[sbox[i] & b] & 0x01;
			}

			sbgen::transform_utils::fwht_transform(truth_table, spectre);

			for (int i = 0; i < 256; i++)
			{
				int spectre_temp = spectre[i];

				if (spectre_temp < 0)
					spectre_temp = -spectre_temp;

				if (spectre_temp <= data->x)
					continue;

				int32_t val = ((spectre_temp - data->y) >= 0) ? 
					(spectre_temp - data->y) : -(spectre_temp - data->y);

				double part = (val>>2)*data->r;
				cost.cost += pow(2.0, part);

				if (spectre_temp > max_spectre)
					max_spectre = static_cast<int32_t>(spectre_temp);
			}
		}

		cost.nonlinearity = 128 - max_spectre / 2;

		return cost;
	}

}; // namespace sbgen

#endif // _COST_FUNCTION_H_
