#ifndef  _SBGEN_TEST_UTILS_H_
#define _SBGEN_TEST_UTILS_H_

#include <cstdio>

#define print_error(...) \
	do {\
		printf("\033[0;31m"); \
		printf("[ERROR] "); \
		printf("\033[0m"); \
		printf(__VA_ARGS__);\
	} while(false);
    
#define print_ok(...) \
	do {\
		printf("\033[0;32m"); \
		printf("[OK] "); \
		printf("\033[0m"); \
		printf(__VA_ARGS__);\
	} while(false);
    
#define print_ignored(...) \
	do {\
		printf("\033[0;33m"); \
		printf("[IGNORED] "); \
		printf("\033[0m"); \
		printf(__VA_ARGS__);\
	} while(false);
	
#define run_test(test_func, test_name, result)\
	do {\
		if(test_func()) {\
			result = result + 1;\
			print_ok(test_name)\
		} else\
			print_error(test_name)\
	} while(false)
	
#define ACCEPT_TEST 1
#define REJECT_TEST  0

#endif
