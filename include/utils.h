#ifndef _SBOX_UTILS_H_
#define _SBOX_UTILS_H_

#include <stdio.h>
#include <stdarg.h>
#include <array>

class logger {
public:
	bool is_enabled = false;
	static logger* me() {
		static logger me;
		return &me;
	}

	void log(const char* fmt, ...) {
		if (!is_enabled)
			return;

		va_list argptr;
		va_start(argptr, fmt);
		vfprintf(stderr, fmt, argptr);
		va_end(argptr);
	}

	void print_sbox(std::array<uint8_t, 256> &sbox) {
		printf("target sbox:\n");
		for (int i = 0;i < 256;i++) {
			if (i != 0 && i % 16 == 0)
				printf("\n");
			printf("0x%02X, ", sbox[i]);
		}
		printf("\n");
	}
};

#define LOG_ON logger::me()->is_enabled = true;
#define LOG_OFF logger::me()->is_enabled = false;
#define LOG(_fmt, ...) logger::me()->log(_fmt, __VA_ARGS__)
#define PRINT_SBOX(sbox) logger::me()->print_sbox(sbox)




#endif