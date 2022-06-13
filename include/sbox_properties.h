// Implementation of properties of s-box and useful transformations
#ifndef  _SBOX_PROPERTIES_H_
#define _SBOX_PROPERTIES_H_

#include <random>

namespace sbgen {

/**
 * @brief Holds S-box transformations
 *
 * Contains set of functions and precomuted values for property calculation
 **/
class transform_utils {
public:
	
  /**
   * @brief Walsh-Hadamard transform of boolean function
   *
   * Compute Walsh-Hadamard transform and return spectre via parameter
   *
   * @param truth_table
   *    thruth table of boolean function
   * @param spectre
   *    output param. Spectre of boolean function
   *@returns
   *   spectre of boolean function via second parameter
   */
	static void fwht_transform
	(
		uint8_t* truth_table, 
		int32_t* spectre
	)
	{
		int32_t step = 1;
		
		for (int32_t i = 0; i < 256; i++)
			spectre[i] = -2 * truth_table[255 - i] + 1;
		
		while (step < 256) 
		{
			int32_t left = 0;
			int32_t numOfBlocks = (256 / (step * 2));
			
			for (int32_t i = 0; i < numOfBlocks; i++) 
			{
				int32_t right = left + step;

				for (int32_t j = 0; j < step; j++) 
				{
					int32_t a = spectre[right];
					int32_t b = spectre[left];
					
					spectre[left] = a + b;
					spectre[right] = a - b;
					
					left++;
					right++;
				}
				left = right;
			}
			step *= 2;
		}
		return;
	}

	
  /**
   * @brief Compute values of all monomials degree less of equal max_deg for set of boolean values
   * 
   * This function used in algebraic immunity calculation.
   * 
   * @param x
   *    initial values
   * @param monomials
   *    monomial values
   * @param max_deg
   *    maximal degree of monomials
   *@returns
   *   calculated monomial values for set x via second parameter
   */
	static void to_monomials(bool* x, bool* monomials, int max_deg) 
	{
		monomials[0] = 1;
		//monomials x1,x8,y1,...,y8
		for (int i = 1;i <= 16;i++)
			monomials[i] = x[i - 1];

		if (max_deg < 2)
			return;
		
		uint32_t pos = 17;
		//monomials x1x2
		for (int i = 1;i < 16;i++)
		{
			for (int j = i + 1;j <= 16;j++) 
			{
				monomials[pos] = monomials[i] & monomials[j];
				pos++;
			}
		}
		if (max_deg < 3)
			return;
    
		//monomials x1x2x3
		for (int i = 1;i < 15;i++) 
		{
			for (int j = i + 1;j <= 16;j++) 
			{
				for (int k = j + 1;k <= 16;k++) 
				{
					monomials[pos] = monomials[i] & monomials[j] & monomials[k];
					pos++;
				}
			}
		}
	}


  /**
   * @brief Compute matrix rank
   * 
   * Compute matrix rank via gauss elimination for 137x256 (or less) matrix with elements in GF(2).
   * Used in algebraic degree computation.
   * 
   * @param a
   *    matrix (system of linear equations)
   * @param n
   *    columns count
   * @param m
   *    rows count
   *@returns
   *   rank of matix a
   */
	static uint32_t gauss_elimination
	(
		bool a[137][256],
		int n,
		int m
	) 
	{
		uint32_t rank = 256;
		bool line_used[137];
		
		for (int i = 0;i < 137;i++)
			line_used[i] = false;
		
		for (int i = 0; i < m; ++i) 
		{
			int j;
			
			for (j = 0; j < n; ++j)
				if (!line_used[j] && a[j][i])
					break;
				
			if (j == n)
				--rank;
			
			else 
			{
				line_used[j] = true;
				
				for (int k = 0; k < n; ++k)
					if (k != j && a[k][i])
						for (int p = i + 1; p < m; ++p)
							a[k][p] = !a[k][p] != !(a[j][p] && a[k][i]);//logical xor: !a != !b
			}
		}
		return rank;
	}
	
	/**
	* @brief Erase fixed and inverse fixed points in sbox 
	* 
	* @param sbox
	*    matrix (system of linear equations)
	* @param seed
	*/
	static void erase_fixed_points(std::array<uint8_t, 256>& sbox, int32_t seed=0xdeadbeef) 
	{
		std::mt19937 gen(seed);
		std::uniform_int_distribution<int> distrib(0, 255);
		std::array<uint8_t, 256> sbox_new;
		uint8_t shift1;
		uint8_t shift2;
		uint8_t mask;
		uint8_t flag;
	
		while (1) {
		
			shift1 = distrib(gen) % 7 + 1;
			shift2 = distrib(gen) % 7 + 1;

			mask = distrib(gen);

			for (uint32_t i = 0;i < 256;i++) 
			{
				uint8_t new_i = (((i << shift1) & 0xFF) | ((i >> (8 - shift1))));
				sbox_new[i] = (((sbox[new_i] << shift2) & 0xFF) | ((sbox[new_i] >> (8 - shift2))));
				sbox_new[i] ^= mask;
			}
		
			flag = 0;
			for (uint32_t i = 0;i < 256;i++) 
			{
				sbox[i] = sbox_new[i];
				if (sbox[i] == i || sbox[i] == (i ^ 0xFF)) 
				{
					flag = 1;
				}
			}
		
			if (flag == 0)
				break;
		}
		return;
	}
	
	// count of "1" in binary representation on numbers 0 - 255
	static  constexpr uint8_t one_bits[256] = {
		0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 
		3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 
		3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 
		4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 
		3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 
		6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 
		4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 
		6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
	};
	
	// bit representation of 0-255
	static constexpr uint8_t bits[256][8] = {
		{0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,1}, {0,0,0,0,0,0,1,0}, {0,0,0,0,0,0,1,1},
		{0,0,0,0,0,1,0,0}, {0,0,0,0,0,1,0,1}, {0,0,0,0,0,1,1,0}, {0,0,0,0,0,1,1,1},
		{0,0,0,0,1,0,0,0}, {0,0,0,0,1,0,0,1}, {0,0,0,0,1,0,1,0}, {0,0,0,0,1,0,1,1},
		{0,0,0,0,1,1,0,0}, {0,0,0,0,1,1,0,1}, {0,0,0,0,1,1,1,0}, {0,0,0,0,1,1,1,1},
		{0,0,0,1,0,0,0,0}, {0,0,0,1,0,0,0,1}, {0,0,0,1,0,0,1,0}, {0,0,0,1,0,0,1,1},
		{0,0,0,1,0,1,0,0}, {0,0,0,1,0,1,0,1}, {0,0,0,1,0,1,1,0}, {0,0,0,1,0,1,1,1},
		{0,0,0,1,1,0,0,0}, {0,0,0,1,1,0,0,1}, {0,0,0,1,1,0,1,0}, {0,0,0,1,1,0,1,1},
		{0,0,0,1,1,1,0,0}, {0,0,0,1,1,1,0,1}, {0,0,0,1,1,1,1,0}, {0,0,0,1,1,1,1,1},
		{0,0,1,0,0,0,0,0}, {0,0,1,0,0,0,0,1}, {0,0,1,0,0,0,1,0}, {0,0,1,0,0,0,1,1},
		{0,0,1,0,0,1,0,0}, {0,0,1,0,0,1,0,1}, {0,0,1,0,0,1,1,0}, {0,0,1,0,0,1,1,1},
		{0,0,1,0,1,0,0,0}, {0,0,1,0,1,0,0,1}, {0,0,1,0,1,0,1,0}, {0,0,1,0,1,0,1,1},
		{0,0,1,0,1,1,0,0}, {0,0,1,0,1,1,0,1}, {0,0,1,0,1,1,1,0}, {0,0,1,0,1,1,1,1},
		{0,0,1,1,0,0,0,0}, {0,0,1,1,0,0,0,1}, {0,0,1,1,0,0,1,0}, {0,0,1,1,0,0,1,1},
		{0,0,1,1,0,1,0,0}, {0,0,1,1,0,1,0,1}, {0,0,1,1,0,1,1,0}, {0,0,1,1,0,1,1,1},
		{0,0,1,1,1,0,0,0}, {0,0,1,1,1,0,0,1}, {0,0,1,1,1,0,1,0}, {0,0,1,1,1,0,1,1},
		{0,0,1,1,1,1,0,0}, {0,0,1,1,1,1,0,1}, {0,0,1,1,1,1,1,0}, {0,0,1,1,1,1,1,1},
		{0,1,0,0,0,0,0,0}, {0,1,0,0,0,0,0,1}, {0,1,0,0,0,0,1,0}, {0,1,0,0,0,0,1,1},
		{0,1,0,0,0,1,0,0}, {0,1,0,0,0,1,0,1}, {0,1,0,0,0,1,1,0}, {0,1,0,0,0,1,1,1},
		{0,1,0,0,1,0,0,0}, {0,1,0,0,1,0,0,1}, {0,1,0,0,1,0,1,0}, {0,1,0,0,1,0,1,1},
		{0,1,0,0,1,1,0,0}, {0,1,0,0,1,1,0,1}, {0,1,0,0,1,1,1,0}, {0,1,0,0,1,1,1,1},
		{0,1,0,1,0,0,0,0}, {0,1,0,1,0,0,0,1}, {0,1,0,1,0,0,1,0}, {0,1,0,1,0,0,1,1},
		{0,1,0,1,0,1,0,0}, {0,1,0,1,0,1,0,1}, {0,1,0,1,0,1,1,0}, {0,1,0,1,0,1,1,1},
		{0,1,0,1,1,0,0,0}, {0,1,0,1,1,0,0,1}, {0,1,0,1,1,0,1,0}, {0,1,0,1,1,0,1,1},
		{0,1,0,1,1,1,0,0}, {0,1,0,1,1,1,0,1}, {0,1,0,1,1,1,1,0}, {0,1,0,1,1,1,1,1},
		{0,1,1,0,0,0,0,0}, {0,1,1,0,0,0,0,1}, {0,1,1,0,0,0,1,0}, {0,1,1,0,0,0,1,1},
		{0,1,1,0,0,1,0,0}, {0,1,1,0,0,1,0,1}, {0,1,1,0,0,1,1,0}, {0,1,1,0,0,1,1,1},
		{0,1,1,0,1,0,0,0}, {0,1,1,0,1,0,0,1}, {0,1,1,0,1,0,1,0}, {0,1,1,0,1,0,1,1},
		{0,1,1,0,1,1,0,0}, {0,1,1,0,1,1,0,1}, {0,1,1,0,1,1,1,0}, {0,1,1,0,1,1,1,1},
		{0,1,1,1,0,0,0,0}, {0,1,1,1,0,0,0,1}, {0,1,1,1,0,0,1,0}, {0,1,1,1,0,0,1,1},
		{0,1,1,1,0,1,0,0}, {0,1,1,1,0,1,0,1}, {0,1,1,1,0,1,1,0}, {0,1,1,1,0,1,1,1},
		{0,1,1,1,1,0,0,0}, {0,1,1,1,1,0,0,1}, {0,1,1,1,1,0,1,0}, {0,1,1,1,1,0,1,1},
		{0,1,1,1,1,1,0,0}, {0,1,1,1,1,1,0,1}, {0,1,1,1,1,1,1,0}, {0,1,1,1,1,1,1,1},
		{1,0,0,0,0,0,0,0}, {1,0,0,0,0,0,0,1}, {1,0,0,0,0,0,1,0}, {1,0,0,0,0,0,1,1},
		{1,0,0,0,0,1,0,0}, {1,0,0,0,0,1,0,1}, {1,0,0,0,0,1,1,0}, {1,0,0,0,0,1,1,1},
		{1,0,0,0,1,0,0,0}, {1,0,0,0,1,0,0,1}, {1,0,0,0,1,0,1,0}, {1,0,0,0,1,0,1,1},
		{1,0,0,0,1,1,0,0}, {1,0,0,0,1,1,0,1}, {1,0,0,0,1,1,1,0}, {1,0,0,0,1,1,1,1},
		{1,0,0,1,0,0,0,0}, {1,0,0,1,0,0,0,1}, {1,0,0,1,0,0,1,0}, {1,0,0,1,0,0,1,1},
		{1,0,0,1,0,1,0,0}, {1,0,0,1,0,1,0,1}, {1,0,0,1,0,1,1,0}, {1,0,0,1,0,1,1,1},
		{1,0,0,1,1,0,0,0}, {1,0,0,1,1,0,0,1}, {1,0,0,1,1,0,1,0}, {1,0,0,1,1,0,1,1},
		{1,0,0,1,1,1,0,0}, {1,0,0,1,1,1,0,1}, {1,0,0,1,1,1,1,0}, {1,0,0,1,1,1,1,1},
		{1,0,1,0,0,0,0,0}, {1,0,1,0,0,0,0,1}, {1,0,1,0,0,0,1,0}, {1,0,1,0,0,0,1,1},
		{1,0,1,0,0,1,0,0}, {1,0,1,0,0,1,0,1}, {1,0,1,0,0,1,1,0}, {1,0,1,0,0,1,1,1},
		{1,0,1,0,1,0,0,0}, {1,0,1,0,1,0,0,1}, {1,0,1,0,1,0,1,0}, {1,0,1,0,1,0,1,1},
		{1,0,1,0,1,1,0,0}, {1,0,1,0,1,1,0,1}, {1,0,1,0,1,1,1,0}, {1,0,1,0,1,1,1,1},
		{1,0,1,1,0,0,0,0}, {1,0,1,1,0,0,0,1}, {1,0,1,1,0,0,1,0}, {1,0,1,1,0,0,1,1},
		{1,0,1,1,0,1,0,0}, {1,0,1,1,0,1,0,1}, {1,0,1,1,0,1,1,0}, {1,0,1,1,0,1,1,1},
		{1,0,1,1,1,0,0,0}, {1,0,1,1,1,0,0,1}, {1,0,1,1,1,0,1,0}, {1,0,1,1,1,0,1,1},
		{1,0,1,1,1,1,0,0}, {1,0,1,1,1,1,0,1}, {1,0,1,1,1,1,1,0}, {1,0,1,1,1,1,1,1},
		{1,1,0,0,0,0,0,0}, {1,1,0,0,0,0,0,1}, {1,1,0,0,0,0,1,0}, {1,1,0,0,0,0,1,1},
		{1,1,0,0,0,1,0,0}, {1,1,0,0,0,1,0,1}, {1,1,0,0,0,1,1,0}, {1,1,0,0,0,1,1,1},
		{1,1,0,0,1,0,0,0}, {1,1,0,0,1,0,0,1}, {1,1,0,0,1,0,1,0}, {1,1,0,0,1,0,1,1},
		{1,1,0,0,1,1,0,0}, {1,1,0,0,1,1,0,1}, {1,1,0,0,1,1,1,0}, {1,1,0,0,1,1,1,1},
		{1,1,0,1,0,0,0,0}, {1,1,0,1,0,0,0,1}, {1,1,0,1,0,0,1,0}, {1,1,0,1,0,0,1,1},
		{1,1,0,1,0,1,0,0}, {1,1,0,1,0,1,0,1}, {1,1,0,1,0,1,1,0}, {1,1,0,1,0,1,1,1},
		{1,1,0,1,1,0,0,0}, {1,1,0,1,1,0,0,1}, {1,1,0,1,1,0,1,0}, {1,1,0,1,1,0,1,1},
		{1,1,0,1,1,1,0,0}, {1,1,0,1,1,1,0,1}, {1,1,0,1,1,1,1,0}, {1,1,0,1,1,1,1,1},
		{1,1,1,0,0,0,0,0}, {1,1,1,0,0,0,0,1}, {1,1,1,0,0,0,1,0}, {1,1,1,0,0,0,1,1},
		{1,1,1,0,0,1,0,0}, {1,1,1,0,0,1,0,1}, {1,1,1,0,0,1,1,0}, {1,1,1,0,0,1,1,1},
		{1,1,1,0,1,0,0,0}, {1,1,1,0,1,0,0,1}, {1,1,1,0,1,0,1,0}, {1,1,1,0,1,0,1,1},
		{1,1,1,0,1,1,0,0}, {1,1,1,0,1,1,0,1}, {1,1,1,0,1,1,1,0}, {1,1,1,0,1,1,1,1},
		{1,1,1,1,0,0,0,0}, {1,1,1,1,0,0,0,1}, {1,1,1,1,0,0,1,0}, {1,1,1,1,0,0,1,1},
		{1,1,1,1,0,1,0,0}, {1,1,1,1,0,1,0,1}, {1,1,1,1,0,1,1,0}, {1,1,1,1,0,1,1,1},
		{1,1,1,1,1,0,0,0}, {1,1,1,1,1,0,0,1}, {1,1,1,1,1,0,1,0}, {1,1,1,1,1,0,1,1},
		{1,1,1,1,1,1,0,0}, {1,1,1,1,1,1,0,1}, {1,1,1,1,1,1,1,0}, {1,1,1,1,1,1,1,1}
	};
}; // class transform_utils 

// Macro definitions for s-box properties
#define SBGEN_MAX_PROPERTIES_NUMBER	3

#define SBGEN_NONLINEARITY				0
#define SBGEN_DELTA_UNIFORMITY		1
#define SBGEN_ALGEBRAIC_IMMUNITY	2

#define SBGEN_USE_NONLINEARITY_FLAG					(1<<0)
#define SBGEN_USE_DELTA_UNIFORMITY_FLAG			(1<<1)
#define SBGEN_USE_ALGEBRAIC_IMMUNITY_FLAG 	(1<<2)

#define setup_property(config, property_index, target_value)\
	do {\
		(config)->properties_config |= (1<<property_index);\
		(config)->target_properties[property_index] = target_value;\
	} while(false);

/**
 * @brief S-box properties
 *
 * Contains set of functions that compute s-box properties
 **/
class properties {
public:

	/**
	* @brief Is fixed points or inverse fixed points in sbox
	* 
	* @param sbox
	*    target sbox
	*@returns
	*   true - fixed points present, otherwise false
	*/
	static int32_t fixed_points(std::array<uint8_t, 256>& sbox)
	{
		int32_t flag = 0;
		for (uint32_t i = 0;i < 256;i++) 
		{
			if (sbox[i] == i || sbox[i] == (i ^ 0xFF)) 
			{
				flag = 1;
				break;
			}
		}
		
		return flag;
	}
	
	
	/**
	* @brief Nonlinearity
	* 
	* Compute S-box nonlinearity
	* 
	* @param sbox
	*    target sbox
	*@returns
	*   sbox nonlinearity
	*/
	static int32_t nonlinearity(std::array<uint8_t, 256>& sbox)
	{
		uint8_t     truth_table[256];	
		int         spectre[256];
		int         max_spectre = 0;
    
		max_spectre = 0;
		
		for (uint16_t b = 1; b <256 ; b++) 
		{
			for (int i = 0;i < 256;i++) {
				truth_table[i] = transform_utils::one_bits[sbox[i] & b] & 0x01;
			}

			transform_utils::fwht_transform(truth_table, spectre);

			for (int i = 0;i < 256;i++) 
			{
				if (spectre[i] < 0)
					spectre[i] = -spectre[i];

				if (spectre[i] > max_spectre)
					max_spectre = spectre[i];
			}
		}

		return 128 - max_spectre / 2;
	}
	
  /**
   * @brief Delta uniformity
   * 
   * Compute S-box delta uniformity
   * 
   * @param sbox
   *    target sbox
   *@returns
   *   sbox delta uniformity
   */
	static int32_t  delta_uniformity(std::array<uint8_t, 256>& sbox) 
	{
		int32_t max_res = 0;
		int32_t res = 0;

		for (uint32_t a = 1;a < 256;a++) 
		{
			for (uint32_t b = 0;b < 256;b++) 
			{
				res = 0;
				for (uint32_t x = 0;x < 256;x++)
					if ((sbox[x] ^ sbox[x ^ a]) == b)
						res++;
                
				if (res > max_res)
					max_res = res;
			}
		}

		return max_res;
	}
	
  /**
   * @brief Algebraic immunity
   * 
   * Compute S-box algebraic immunity
   * 
   * @param sbox
   *    target sbox
   *@returns
   *   sbox algebraic immunity
   */
	static int32_t algebraic_immunity(std::array<uint8_t, 256>& sbox) 
	{

		
		bool mat[137][256];
		bool tmp[137];
		bool values[16];
		
		for (int i = 0;i < 256;i++) 
		{
			uint8_t y = sbox[i];
			values[0] = transform_utils::bits[i][0]; values[1] = transform_utils::bits[i][1];
			values[2] = transform_utils::bits[i][2]; values[3] = transform_utils::bits[i][3];
			values[4] = transform_utils::bits[i][4]; values[5] = transform_utils::bits[i][5];
			values[6] = transform_utils::bits[i][6]; values[7] = transform_utils::bits[i][7];
			values[8] = transform_utils::bits[y][0]; values[9] = transform_utils::bits[y][1];
			values[10] = transform_utils::bits[y][2]; values[11] = transform_utils::bits[y][3];
			values[12] = transform_utils::bits[y][4]; values[13] = transform_utils::bits[y][5];
			values[14] = transform_utils::bits[y][6]; values[15] = transform_utils::bits[y][7];
			transform_utils::to_monomials((bool*)&values, (bool*)tmp,2);
			for (int j = 0;j < 137;j++)
				mat[j][i] = tmp[j];
		}

		uint32_t rank = transform_utils::gauss_elimination(mat,137,256);

		if (rank == 137)
			return 3;
    
		for (int i = 0;i < 256;i++)
		{
			uint8_t y = sbox[i];
			values[0] = transform_utils::bits[i][0]; values[1] = transform_utils::bits[i][1];
			values[2] = transform_utils::bits[i][2]; values[3] = transform_utils::bits[i][3];
			values[4] = transform_utils::bits[i][4]; values[5] = transform_utils::bits[i][5];
			values[6] = transform_utils::bits[i][6]; values[7] = transform_utils::bits[i][7];
			values[8] = transform_utils::bits[y][0]; values[9] = transform_utils::bits[y][1];
			values[10] = transform_utils::bits[y][2]; values[11] = transform_utils::bits[y][3];
			values[12] = transform_utils::bits[y][4]; values[13] = transform_utils::bits[y][5];
			values[14] = transform_utils::bits[y][6]; values[15] = transform_utils::bits[y][7];
			transform_utils::to_monomials((bool*)&values, (bool*)tmp, 1);
			for (int j = 0;j < 17;j++)
				mat[j][i] = tmp[j];
		}
    
		rank = transform_utils::gauss_elimination(mat, 17, 256);

		if (rank == 17)
			return 2;

		return 1;
	}
	
}; // class properties

    
}; // namespace sbgen

#endif // _SBOX_PROPERTIES_H_
