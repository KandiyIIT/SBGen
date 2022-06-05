# SBGen: 8-bit S-Box generation framework
Header-only multithread c++17 framework for s-box generation. Usage example:

```cpp
    #include "hill_climbing.h"
    #include "utils.h"

	sbgen::hill_climbing_info_t<double> info;

	info.thread_count = 2;
	info.is_log_enabled = false;

	info.try_per_thread = 1000000;
	info.max_frozen_count = 100000;
	
	setup_property( &info, SBGEN_NONLINEARITY, 102);

	
	info.use_random_seed = false;
	info.seed = 0xdeadbeef;



	info.cost_function = sbgen::whs<double>;
	info.cost_data.reset(new sbgen::whs_function_data_t(12, 0));

	//s-box stored in table format in std::optional<std::array<uint8_t,256>>
	auto sbox = sbox::hill_climbing<double>(info);

	if (sbox.has_value())
		PRINT_SBOX(sbox.value());
```
