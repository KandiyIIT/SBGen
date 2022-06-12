// Utils for debug
#ifndef _SBOX_UTILS_H_
#define _SBOX_UTILS_H_

#include <cstdio>
#include <stdarg.h>
#include <array>


/**
 * @brief Simple logging class for debug 
 * 
 * Implemented as a singleton.
 */
class logger {
public:
	bool is_enabled = false;
    
    
  /**
   * @brief Return pointer to actual logger object 
   * 
   * @returns
   *    pointer to actual logger object
   */
	static logger* me() {
		static logger me;
		return &me;
	}

	
  /**
   * @brief Log info to standart output
   *
   * Works just like printf
   *
   * @param fmt
   *    format string
   * 
   * @param ...
   *    list of values for printf
   */
	void log(const char* fmt, ...) {
		if (!is_enabled)
			return;

		va_list argptr;
		va_start(argptr, fmt);
		vfprintf(stderr, fmt, argptr);
		va_end(argptr);
	}

	
  /**
   * @brief Print s-box to standart output
   * 
   * Print s-box in table format 16x16
   *
   * @param sbox
   *    S-box
   */
	void print_sbox(std::array<uint8_t, 256> &sbox) {
		printf("target sbox:\n");
		for (int i = 0;i < 256;i++) {
			if (i != 0 && i % 16 == 0)
				printf("\n");
			printf("0x%02X, ", sbox[i]);
		}
		printf("\n");
	}
	
}; // class logger 


// Some usefull macro definitions
#define LOG_ON logger::me()->is_enabled = true;
#define LOG_OFF logger::me()->is_enabled = false;
#define LOG(_fmt, ...) logger::me()->log(_fmt, __VA_ARGS__)
#define PRINT_SBOX(sbox) logger::me()->print_sbox(sbox)

#endif // _SBOX_UTILS_H_
