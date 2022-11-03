# SBGen: 8-bit S-Box generation framework
Header-only multithread c++17 framework for s-box generation. the library contains implementations of the following generation methods:

* Hill Climbing [[CS04](#CS04)]
* Simulated Annealing [[KRPKK22](#KRPKK22)]

Supported cost functions:

* WHS  [[CS04](#CS04)]
* PCF [[PCR16](#PCR16)]
* WCF [[FAM20](#FAM20)]

The main property for which optimization is performed is non-linearity, but the library also supports a number of other properties:

* Nonlinearity (or Linear approximation probability)
* Delta Uniformity (or Differential approximation probability)
* Algebraic Immunity
* Strict Avalanche Criterion
* Output Bit Independence Criterion

Usage example:

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
	auto sbox = sbgen::hill_climbing<double>(info);

	if (sbox.has_value())
		PRINT_SBOX(sbox.value());
```

See SBGenTool and tests as a source of examples.

# Installation #

The library does not require installation as it contains only header files. We also provide a utility SBGenTool. To compile it, just run the script:

```bash

cd ./SBGenTool
make

```

Interface of SBGenTool:

```
Usage: sbgen --method [METHOD] [OPTIONS] 
List of options:

  --visibility
       Enable verbose mode
  --version
       Print version info
  --help
       Print help message
  --seed
       seed for randomness. Warning: in multithread mode there is
       additional randomnes caused by concurrency
  --method [hill_climbing|simulated_annealing]
       hill_climbing = hill climbing method
       simulated_annealing = simulated annealing method
  --cost_function [whs|wcf|pcf]
       whs = WHS cost function
       wcf = WCF cost function
       pcf = PCF cost function
  --thread_count
       max thread count
  --cost_type [int64_t|double]
       type of variable, where stored s-box cost. Default value - double
  --try_per_thread
       maximal iterations count in method
  --max_frozen_loops
       max iterations count without any chages

Method parameter list:

  --method_params
       params of method in format --method_params={param1,param2,...,paramN}
  hill_climbing:
       Has no free options
  simulated_annealing
       param1: max_outer_loops - maximal outer loop count
       param2: max_inner_loops - maximal inner loop count
       param3: initial_temperature -  initial temperature
       param4: alpha_parameter -  alpha_parameter
       Example: --method_params="{10, 10000, 1000, 0.99}"

Cost function parameter list:

  --cost_function_params
       params of cost function in format --cost_function_params={param1,param2,...,paramN}
  whs
       param1: r
       param2: x
       Example: --cost_function_params="{12, 0}"
  pcf
       param1: n
       Example: --cost_function_params="{5}"
  wcf
       Has no free options

Target properties:

  --nonlinearity
       target nonlinearity value.
  --delta_uniformity
       target delta uniformity value.
  --algebraic_immunity
       target algebraic immunity value.
```

Usage example:

```
sbgen --method=hill_climbing --nonlinearity=104 --try_per_thread=10000000 --thread_count=8 --cost_function=wcf   --delta_uniformity=8 --algebraic_immunity=3 --erase_points --seed=1234 
```
Output:

```
0x65, 0x96, 0x50, 0x66, 0x2E, 0x14, 0x49, 0x23, 0xCC, 0xBC, 0xAE, 0xEE, 0x63, 0xEA, 0x78, 0xD8, 
0x02, 0x05, 0x82, 0x1D, 0xED, 0x9C, 0xF7, 0x33, 0x79, 0xEC, 0xBD, 0x6B, 0x80, 0x9F, 0x97, 0x35, 
0x59, 0x8B, 0xFD, 0x71, 0x20, 0xFF, 0xEB, 0xE2, 0x34, 0x90, 0x28, 0xFC, 0x17, 0x5E, 0xF2, 0xDE, 
0x39, 0x53, 0xB5, 0x0D, 0x84, 0xE8, 0xA1, 0x19, 0x64, 0x4D, 0xFB, 0xCA, 0xB9, 0x56, 0xB2, 0x58, 
0x5F, 0xB8, 0x88, 0xA7, 0xF1, 0x4F, 0x1B, 0x57, 0x0B, 0x74, 0xDC, 0x8D, 0x31, 0xF4, 0x1A, 0x5A, 
0xC6, 0xA4, 0x7F, 0x13, 0x62, 0x75, 0x41, 0xC8, 0xC9, 0x46, 0xE4, 0x07, 0xDB, 0x98, 0x08, 0xB7, 
0x16, 0xF0, 0xAC, 0x67, 0xAA, 0x18, 0x25, 0xD2, 0x60, 0xC3, 0x47, 0x42, 0xAB, 0x1E, 0xA3, 0x6A, 
0x7C, 0x00, 0x36, 0xEF, 0xBB, 0x32, 0x2A, 0x4E, 0xE0, 0x3B, 0x4C, 0x9E, 0x87, 0xF3, 0x95, 0x15, 
0xAF, 0xBA, 0x51, 0xD6, 0xA6, 0x1C, 0x81, 0x61, 0x45, 0x04, 0xE7, 0x55, 0xC1, 0xC0, 0xAD, 0x38, 
0xF9, 0x26, 0xBF, 0x86, 0x89, 0x27, 0x91, 0x03, 0xCE, 0xE1, 0x3A, 0x3D, 0xCB, 0x9A, 0x93, 0xB1, 
0x2C, 0xDA, 0x3E, 0x4A, 0xF6, 0xE9, 0x9D, 0xB4, 0xC2, 0x11, 0xD0, 0x0F, 0x40, 0x4B, 0x6E, 0xA5, 
0xDF, 0x0A, 0x6D, 0x83, 0x8E, 0xE5, 0xC4, 0x7E, 0x54, 0xD9, 0x5B, 0x01, 0x52, 0x7A, 0x43, 0x30, 
0x37, 0x9B, 0xD7, 0x77, 0xD1, 0x92, 0x85, 0xF8, 0x2F, 0xE6, 0xBE, 0x7D, 0xA0, 0xC7, 0x06, 0xFE, 
0x94, 0x8C, 0x3C, 0xB3, 0x10, 0x21, 0x3F, 0x09, 0xB0, 0xC5, 0x99, 0xB6, 0xD4, 0x48, 0x2D, 0x22, 
0xCD, 0x8F, 0x2B, 0x5D, 0xFA, 0x5C, 0x0C, 0x7B, 0x6C, 0xD3, 0x76, 0x0E, 0xF5, 0xA8, 0x6F, 0xD5, 
0x70, 0x29, 0xA9, 0x8A, 0x24, 0x1F, 0xA2, 0x72, 0x68, 0xDD, 0x12, 0x44, 0x69, 0x73, 0xE3, 0xCF, 
NL= 104
DU= 8
AI= 3
Fixed Points= 0
```

<a name="CS04">[CS04]<a/> John A. Clark. The Design of S-Boxes by Simulated Annealing. New Generation Computing 23(3):219-231 (2004)
<a name="KRPKK22">[KRPKK22]<a/> A. Kuznetsov. Optimizing Hill Climbing Algorithm Parameters forGeneration of Cryptographically Strong S-Boxes. https://doi.org/10.21203/rs.3.rs-1657863/v1 (2022)
<a name="PCR16">[PCR16]<a/> S.  Picek,  M.  Cupic,  L.  Rotim.  A  New  Cost  Function  for  Evolution  of  S-Boxes,  Evolutionary Computation. 24 (2016) 695–718. https://doi.org/10.1162/EVCO_a_00191
<a name="FAM20">[FAM20]<a/>  A. Freyre-Echevarría, A. Alanezi, I. Martínez-Díaz, M. Ahmad, A.A. Abd El-Latif, H. Kolivand, A.  Razaq,  An  External  Parameter  Independent  Novel  Cost  Function  for  Evolving  Bijective Substitution-Boxes, Symmetry. 12 (2020) 1896. https://doi.org/10.3390/sym12111896.
