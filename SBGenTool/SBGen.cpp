//=============================================================================

#include <iostream>
#include <functional>
#include <string>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <map>
#include <sstream> 
#include <fstream>

#include "hill_climbing.h"
#include "simulated_annealing.h"
#include "genetic.h"
#include "genetic_internals.h"
#include "cost_function.h"
#include "utils.h"
#include "sbgen_info.h"

//=============================================================================

#define ABORT_MSG(y)                         \
{                                            \
	std::cerr << "sbgen: " << y << std::endl;\
	exit(1);                                 \
}

#define no_argument						0
#define required_argument				1 
#define optional_argument				2

#define version_flag					0x00
#define help_flag						0x01
#define method_flag						0x02
#define cost_function_flag 				0x03
#define thread_count_flag				0x04
#define cost_type_flag					0x05
#define visibility_flag					0x06
#define try_per_thread_flag				0x07
#define max_frozen_loops_flag			0x08
#define method_params_flag				0x09
#define cost_function_params_flag		0x0a
#define nonlinearity_flag				0x0b
#define delta_uniformity_flag			0x0c
#define algebraic_immunity_flag			0x0d
#define seed_flag						0x0e
#define erase_points_flag				0x0f
#define selection_method_flag			0x10
#define crossover_method_flag			0x11
#define sbox_count_flag					0x12
#define to_file_flag					0x13

#define default_try_per_thread			10000
#define default_thread_count			1
#define default_log						0
#define default_erase_points			0
#define default_sbox_count				1

const struct option longopts[] =
{
	{"version",				no_argument,		0,	version_flag},
	{"visibility",			no_argument,		0,	visibility_flag},
	{"help",				no_argument,		0,	help_flag},
	{"erase_points",		no_argument,		0,	erase_points_flag},
	{"method",				required_argument,	0,	method_flag},
	{"cost_function",		required_argument,	0,	cost_function_flag},
	{"thread_count",		required_argument,	0,	thread_count_flag},
	{"cost_type",			required_argument,	0,	cost_type_flag},
	{"try_per_thread",		required_argument,	0,	try_per_thread_flag},
	{"method_params",		required_argument,	0,	method_params_flag},
	{"cost_function_params",required_argument,	0,	cost_function_params_flag},
	{"nonlinearity",		required_argument,	0,	nonlinearity_flag},
	{"delta_uniformity",	required_argument,	0,	delta_uniformity_flag},
	{"algebraic_immunity",	required_argument,	0,	algebraic_immunity_flag},
	{"seed",				required_argument,	0,	seed_flag},
	{"max_frozen_loops",	required_argument,	0,	max_frozen_loops_flag},
	{"selection_method",	required_argument,	0,	selection_method_flag},
	{"crossover_method",	required_argument,	0,	crossover_method_flag},
	{"sbox_count",			required_argument,	0,	sbox_count_flag},
	{"to_file",				required_argument,	0,	to_file_flag},
	{0,0,0,0},
};

using options_t = std::map<int, std::string>;

template<typename T>
using full_cost_info_t =
	std::pair<sbgen::cost_function_t<T>, sbgen::cost_function_data_t*>;

//=============================================================================

void print_help()
{
	std::cout << "Usage: sbgen --method [METHOD] [OPTIONS]\n"
		<< "List of options:\n\n"
		<< "\t--visibility\n"
		<< "\t\tEnable verbose mode\n"
		<< "\t--version\n"
		<< "\tPrint version info\n"
		<< "\t--help\n"
		<< "\t\tPrint help message\n"
		<< "\t--to_file\n"
		<< "\t\tRedirect output to file. (Example: --to_file=\"log.txt\"\n"
		<< "\t--seed\n"
		<< "\t\tseed for randomness. Warning: in multithread mode there is\n"
		<< "\t\tadditional randomnes caused by concurrency\n"
		<< "\t--sbox_count\n"
		<< "\t\ttarget sbox count. Default value - 1.\n"
		<< "\t--method [hill_climbing|simulated_annealing|genetic]\n"
		<< "\t\thill_climbing = hill climbing method\n"
		<< "\t\tsimulated_annealing = simulated annealing method\n"
		<< "\t\tgenetic = genetic method\n"
		<< "\t--cost_function [max_whs|whs|wcf|pcf|cf1|cf2]\n"
		<< "\t\tmax_whs = maxWHS cost function\n"
		<< "\t\twhs = WHS cost function\n"
		<< "\t\twcf = WCF cost function\n"
		<< "\t\tpcf = PCF cost function\n"
		<< "\t\tcf1 = CF1 cost function\n"
		<< "\t\tcf2 = CF2 cost function\n"
		<< "\t--selection_method [basic|rank|roulette]\n"
		<< "\t\tbasic = select only best s-boxes\n"
		<< "\t\trank = rank selection\n"
		<< "\t\troulette = roulette wheel selection\n"
		<< "\t--crossover_method=\"name, count, child\"\n"
		<< "\t\tname = crossover method name [cycle|pmx]\n"
		<< "\t\tcount = crossover pairs count\n"
		<< "\t\tchild = child per parent\n"
		<< "\t\tExample: --crossover_method=\"pmx, 10, 1\"\n"
		<< "\t--thread_count\n"
		<< "\t\tmax thread count\n"
		<< "\t--cost_type [int64_t|double]\n"
		<< "\t\ttype of variable, where stored s-box cost.\n"
		<< "\t\tDefault value - double\n"
		<< "\t--try_per_thread\n"
		<< "\t\tmaximal iterations count in method\n"
		<< "\t--max_frozen_loops\n"
		<< "\t\tmax iterations count without any chages\n"
		<< "\t\t(not actual for genetic)\n\n"
		<< "Method parameter list:\n\n"
		<< "\t--method_params\n"
		<< "\t\tparams of method in format\n"
		<< "\t\t--method_params=\"param1,param2,...,paramN\"\n"
		<< "\thill_climbing:\n"
		<< "\t\tHas no free options\n"
		<< "\tsimulated_annealing\n"
		<< "\t\tparam1: max_outer_loops - maximal outer loop count\n"
		<< "\t\tparam2: max_inner_loops - maximal inner loop count\n"
		<< "\t\tparam3: initial_temperature -  initial temperature\n"
		<< "\t\tparam4: alpha_parameter -  alpha_parameter\n"
		<< "\t\tExample: --method_params=\"10, 10000, 1000, 0.99\"\n"
		<< "\tgenetic\n"
		<< "\t\tparam1: initial_population_count - initial s-box count\n"
		<< "\t\tparam2: mutants_per_parent - mutants count in thread\n"
		<< "\t\tparam3: selection_count - selected s-box count\n"
		<< "\t\tparam4: use_crossover - shoud use crossover?\n"
		<< "\t\tExample: --method_params=\"100, 20, 100, 1\"\n\n"
		<< "Cost function parameter list:\n\n"
		<< "\t--cost_function_params\n"
		<< "\t\tparams of cost function in format\n"
		<< "\t\t--cost_function_params=\"param1,param2,...,paramN\"\n"
		<< "\tmax_whs\n"
		<< "\t\tparam1: r\n"
		<< "\t\tparam2: x\n"
		<< "\t\tExample: --cost_function_params=\"4, 36\"\n"
		<< "\twhs\n"
		<< "\t\tparam1: r\n"
		<< "\t\tparam2: x\n"
		<< "\t\tExample: --cost_function_params=\"12, 0\"\n"
		<< "\tcf1\n"
		<< "\t\tparam1: r\n"
		<< "\t\tparam2: x\n"
		<< "\t\tparam2: y\n"
		<< "\t\tExample: --cost_function_params=\"12, 32, 0\"\n"
		<< "\tcf2\n"
		<< "\t\tparam1: r\n"
		<< "\t\tparam2: x\n"
		<< "\t\tparam2: y\n"
		<< "\t\tExample: --cost_function_params=\"12, 32, 0\"\n"
		<< "\tpcf\n"
		<< "\t\tparam1: n\n"
		<< "\t\tExample: --cost_function_params=\"5\"\n"
		<< "\twcf\n"
		<< "\t\tHas no free options\n\n"
		<< "Target properties:\n\n"
		<< "\t--nonlinearity\n"
		<< "\t\ttarget nonlinearity value.\n"
		<< "\t--delta_uniformity\n"
		<< "\t\ttarget delta uniformity value.\n"
		<< "\t--algebraic_immunity\n"
		<< "\t\ttarget algebraic immunity value.\n"
		<< "\t--erase_fixed_points\n"
		<< "\t\tdelete fixed points via affine transform\n\n"
		<< "Please refer to https://github.com/KandiyIIT/SBGen/README.md\n"
		<< "for more information.\n";
}

//-----------------------------------------------------------------------------

void parse_options(
	int				argc,
	char			**argv,
	options_t&		options)
{
	int				flag;
	int				option_index;
	const char*		short_options = "";

	while ((flag=getopt_long(
		argc, argv, short_options, longopts, &option_index)) != -1)
	{
		std::string parameter_value = "";

		if (flag == '?')
			ABORT_MSG("Unknown parameter");

		if (longopts[option_index].has_arg == required_argument)
			parameter_value = optarg;

		options[longopts[option_index].val] = parameter_value;
	}
}

//-----------------------------------------------------------------------------

void split(
	const std::string			&s,
	char						delim,
	std::vector<std::string>	&elems)
{
	std::stringstream			ss(s);
	std::string					item;

	while (std::getline(ss, item, delim))
	{
		if (item.length() > 0)
		{
			elems.push_back(item);  
		}
	}
}

//-----------------------------------------------------------------------------

std::vector<std::string> split(
	const std::string			&s,
	char						delim)
{
	std::vector<std::string>	elems;

	split(s, delim, elems);

	return elems;
}

//-----------------------------------------------------------------------------

template<typename T>
full_cost_info_t<T> get_cost_function(
	options_t&			options)
{
	auto parameter = options.find(cost_function_flag);

	if (parameter != options.end())
	{
		if (parameter->second == "whs")
		{
			int32_t x = 0;
			int32_t r = 3;

			parameter = options.find(cost_function_params_flag);

			if (parameter != options.end()) 
			{
				std::vector<std::string> values =
					split(parameter->second, ',');
				
				if (values.size() != 2)
				{
					ABORT_MSG("Invalid parameters count for whs function");
				}

				try {
					r = std::stoi(values[0]);
					x = std::stoi(values[1]);
				}
				catch (const std::invalid_argument & e) 
				{
					ABORT_MSG(e.what());
				}
				catch (const std::out_of_range & e) 
				{
					ABORT_MSG(e.what());
				}
				
				return std::make_pair<sbgen::cost_function_t<T>,
					sbgen::cost_function_data_t*>(sbgen::whs<T>,
					new sbgen::whs_function_data_t(r, x));
			} 
			else
			{
				ABORT_MSG("Can't find whs parameters");
			}
		} else if (parameter->second == "max_whs")
		{
			int32_t x = 0;
			int32_t r = 3;

			parameter = options.find(cost_function_params_flag);

			if (parameter != options.end()) 
			{
				std::vector<std::string> values =
					split(parameter->second, ',');
				
				if (values.size() != 2)
				{
					ABORT_MSG("Invalid parameters count for whs function");
				}

				try {
					r = std::stoi(values[0]);
					x = std::stoi(values[1]);
				}
				catch (const std::invalid_argument & e) 
				{
					ABORT_MSG(e.what());
				}
				catch (const std::out_of_range & e) 
				{
					ABORT_MSG(e.what());
				}
				
				return std::make_pair<sbgen::cost_function_t<T>,
					sbgen::cost_function_data_t*>(sbgen::max_whs<T>,
					new sbgen::max_whs_function_data_t(r, x));
			} 
			else
			{
				ABORT_MSG("Can't find max_whs parameters");
			}
		} else if (parameter->second == "cf1")
		{
			int32_t x = 0;
			int32_t y = 0;
			int32_t r = 0;

			parameter = options.find(cost_function_params_flag);

			if (parameter != options.end()) 
			{
				std::vector<std::string> values =
					split(parameter->second, ',');

				if (values.size() != 3)
				{
					ABORT_MSG("Invalid parameters count for cf1 function");
				}

				try {
					r = std::stoi(values[0]);
					x = std::stoi(values[1]);
					y = std::stoi(values[2]);
				}
				catch (const std::invalid_argument & e) 
				{
					ABORT_MSG(e.what());
				}
				catch (const std::out_of_range & e) 
				{
					ABORT_MSG(e.what());
				}

				return std::make_pair<sbgen::cost_function_t<T>,
					sbgen::cost_function_data_t*>(sbgen::cf1<T>,
					new sbgen::cf1_function_data_t(r, x, y));
			} 
			else
			{
				ABORT_MSG("Can't find cf1 parameters");
			}
		} else if (parameter->second == "cf2")
		{
			int32_t x = 0;
			int32_t y = 0;
			int32_t r = 0;

			parameter = options.find(cost_function_params_flag);

			if (parameter != options.end()) 
			{
				std::vector<std::string> values =
					split(parameter->second, ',');
				
				if(values.size()!=3)
				{
					ABORT_MSG("Invalid parameters count for cf2 function");
				}
				
				try {
					r = std::stoi(values[0]);
					x = std::stoi(values[1]);
					y = std::stoi(values[2]);
				}
				catch (const std::invalid_argument & e) 
				{
					ABORT_MSG(e.what());
				}
				catch (const std::out_of_range & e) 
				{
					ABORT_MSG(e.what());
				}
				return std::make_pair<sbgen::cost_function_t<T>,
					sbgen::cost_function_data_t*>(sbgen::cf2<T>,
					new sbgen::cf2_function_data_t(r, x, y));
			} 
			else
			{
				ABORT_MSG("Can't find cf2 parameters");
			}
		} else if(parameter->second == "pcf")
		{
			int32_t n = 5;

			parameter = options.find(cost_function_params_flag);

			if (parameter != options.end()) 
			{
				std::vector<std::string> values =
					split(parameter->second, ',');

				if(values.size()!=1)
				{
					ABORT_MSG("Invalid parameters count for pcf function");
				}
				
				try {
					n = std::stoi(values[0]);
				}
				catch (const std::invalid_argument & e) 
				{
					ABORT_MSG(e.what());
				}
				catch (const std::out_of_range & e) 
				{
					ABORT_MSG(e.what());
				}
				
				return std::make_pair<sbgen::cost_function_t<T>,
					sbgen::cost_function_data_t*>(sbgen::pcf<T>,
					new sbgen::pcf_function_data_t(n));
			} 
			else
			{
				ABORT_MSG("Can't find pcf parameters");
			}
		}
		else if(parameter->second == "wcf")
		{
			return std::make_pair<sbgen::cost_function_t<T>,
				sbgen::cost_function_data_t*>(sbgen::wcf<T>,
				new sbgen::wcf_function_data_t());
		}

		ABORT_MSG("Unknown cost function. See help for avaliable functions");
	}

	return std::make_pair<sbgen::cost_function_t<T>,
		sbgen::cost_function_data_t*>(
			sbgen::wcf<T>,new sbgen::wcf_function_data_t());
}

//-----------------------------------------------------------------------------

void setup_properties(
	sbgen::properties_info_t*	info,
	options_t&					options) 
{
	info->properties_config = 0;

	try 
	{
		auto property = options.find(nonlinearity_flag);
		if (property != options.end()) 
		{
			setup_property( info,
				SBGEN_NONLINEARITY, std::stoi(property->second));
		}
		else
		{
			ABORT_MSG("Need target nonlinearity");
		}
		
		property = options.find(delta_uniformity_flag);
		if (property != options.end()) 
		{
			setup_property( info,
				SBGEN_DELTA_UNIFORMITY, std::stoi(property->second));
		}
		
		property = options.find(algebraic_immunity_flag);
		if (property != options.end()) 
		{
			setup_property( info,
				SBGEN_ALGEBRAIC_IMMUNITY, std::stoi(property->second));
		}
		
		property = options.find(seed_flag);
		if (property != options.end()) 
		{
				info->use_random_seed = false;
				info->seed = std::stoi(property->second);
		}
		else
		{
				info->use_random_seed = true;
		}
	}
	catch (const std::invalid_argument & e) 
	{
		ABORT_MSG(e.what());
	}
	catch (const std::out_of_range & e) 
	{
		ABORT_MSG(e.what());
	}
}

//-----------------------------------------------------------------------------

template<typename T>
auto get_selection_method(
	sbgen::genetic_info_t<T>&	info,
	options_t&					options)
{
	auto parameter = options.find(selection_method_flag);

	if (parameter != options.end()) 
	{
		if (parameter->second == "basic")
			return sbgen::selectors::basic_selection<T>;

		if (parameter->second == "rank")
			return sbgen::selectors::rank_sequential_selection<T>;

		if (parameter->second == "roulette")
			return sbgen::selectors::roulette_wheel_sequential_selection<T>;

		ABORT_MSG("Unknown selection method. See help.");
	}

	return sbgen::selectors::basic_selection<T>;
}

//-----------------------------------------------------------------------------

template<typename T>
void setup_crossover_properties(
	sbgen::genetic_info_t<T>&	info,
	options_t&					options)
{
	if (info.use_crossover == false)
		return;

	auto parameter = options.find(crossover_method_flag);

	if (parameter != options.end()) 
	{
		std::vector<std::string> values = split(parameter->second, ',');
		if(values.size() != 3)
			ABORT_MSG("Invalid crossover parameters. See help");
		
		try 
		{
			std::string crossover_name = values[0];
			info.crossover_count = std::stoi(values[1]);
			info.child_per_parent = std::stoi(values[2]);

			if (crossover_name == "cycle")
				info.crossover_method = sbgen::cossovers::cycle;

			else if (crossover_name == "pmx")
				info.crossover_method = sbgen::cossovers::pmx;

			else ABORT_MSG("Unknown crossover method. See help.");
		}
		catch (const std::invalid_argument & e) 
		{
			ABORT_MSG(e.what());
		}
		catch (const std::out_of_range & e) 
		{
			ABORT_MSG(e.what());
		}
	}
	else
	{
		ABORT_MSG("Can't find genetic crossover parameters");
	}
}

//-----------------------------------------------------------------------------

template<typename T>
void run_generator(
	options_t& options)
{
	int32_t thread_count = default_thread_count;
	int32_t try_per_thread = default_try_per_thread;
	int32_t max_frozen_count = try_per_thread;
	int32_t visibility = default_log;
	int32_t erase_points = default_erase_points;
	
	auto parameter = options.find(method_flag);
	if (parameter != options.end()) 
	{
		try 
		{
			auto property = options.find(thread_count_flag);
			if (property != options.end())
				thread_count = std::stoi(property->second);
			
			property = options.find(try_per_thread_flag);
			if (property != options.end())
				try_per_thread = std::stoi(property->second);
			
			property = options.find(max_frozen_loops_flag);
			if (property != options.end())
				max_frozen_count = std::stoi(property->second);
			else
				max_frozen_count = try_per_thread;
            
			property = options.find(visibility_flag);
			if (property != options.end())
				visibility = 1;
			
			property = options.find(erase_points_flag);
			if (property != options.end())
				erase_points = 1;

		}
		catch (const std::invalid_argument & e) 
		{
			ABORT_MSG(e.what());
		}
		catch (const std::out_of_range & e) 
		{
			ABORT_MSG(e.what());
		}
		
		if(parameter->second == "hill_climbing")
		{
			auto cost = get_cost_function<T>(options);
			sbgen::hill_climbing_info_t<T> info;
			
			info.thread_count = thread_count;
			info.try_per_thread = try_per_thread;
			info.max_frozen_count = max_frozen_count;
			info.cost_function = cost.first;
			info.cost_data.reset(cost.second);
			info.is_log_enabled = visibility;
			setup_properties(
				static_cast<sbgen::properties_info_t*>(&info), options);
			
			std::cout<<"Starting hill climbing..."<<std::endl;
			std::cout<<"Parameters:"<<std::endl;
			std::cout<<"Thread count: "<<info.thread_count<<std::endl;
			std::cout<<"Try per thread: "<<info.try_per_thread<<std::endl;
			std::cout<<"Max frozen loops: "<<info.max_frozen_count<<std::endl;
			std::cout<<"Log level: "<<info.is_log_enabled<<std::endl;
			std::cout<<"Cost Function: "<<info.cost_data->name()<<std::endl;

			if(info.properties_config & SBGEN_USE_NONLINEARITY_FLAG)
			{
				std::cout<<"target NL: "
					<<info.target_properties[SBGEN_NONLINEARITY]
					<<std::endl;
			} else std::cout<<"NL not used\n";

			if(info.properties_config & SBGEN_USE_DELTA_UNIFORMITY_FLAG)
			{
				std::cout<<"target DU: "
					<<info.target_properties[SBGEN_DELTA_UNIFORMITY]
					<<std::endl;
			}else std::cout<<"DU not used\n";

			if(info.properties_config & SBGEN_USE_ALGEBRAIC_IMMUNITY_FLAG)
			{
				std::cout<<"target AI: "
					<<info.target_properties[SBGEN_ALGEBRAIC_IMMUNITY]
					<<std::endl;
			}else std::cout<<"AI not used\n";

			if(info.use_random_seed==true)
				std::cout<<"Seed: random"<<std::endl;
			else
				std::cout<<"Seed: "<<info.seed<<std::endl;

			auto sbox = sbgen::hill_climbing<T>(info);

			if (sbox.has_value()) 
			{
				auto sb = sbox.value();
				
				if(erase_points) 
				{
					sbgen::transform_utils::erase_fixed_points(sb, info.seed);
				}

				PRINT_SBOX(sb);
				
				std::cout<<"NL= "<<
					sbgen::properties::nonlinearity(sb)<<std::endl;
				std::cout<<"DU= "<<
					sbgen::properties::delta_uniformity(sb)<<std::endl;
				std::cout<<"AI= "
					<<sbgen::properties::algebraic_immunity(sb)<<std::endl;
				std::cout<<"Fixed Points= "
					<<sbgen::properties::fixed_points(sb)<<std::endl;
			} else {
				ABORT_MSG("SBox not found. Try another parameters");
			}
		} 
		else if(parameter->second == "simulated_annealing")
		{
			auto cost = get_cost_function<T>(options);
			sbgen::simulated_annealing_info_t<T> info;
			info.thread_count = thread_count;
			info.try_per_thread = try_per_thread;
			info.max_frozen_outer_loops = max_frozen_count;
			info.cost_function = cost.first;
			info.cost_data.reset(cost.second);
			info.is_log_enabled = visibility;
			setup_properties(
				static_cast<sbgen::properties_info_t*>(&info), options);

			parameter = options.find(method_params_flag);
			if (parameter != options.end()) 
			{
				std::vector<std::string> values =
					split(parameter->second, ',');
				if(values.size() != 4)
					ABORT_MSG("Invalid simulated annealing parameters.");

				try 
				{
					info.max_outer_loops = std::stoi(values[0]);
					info.max_inner_loops = std::stoi(values[1]);
					info.initial_temperature = std::stod(values[2]);
					info.alpha_parameter = std::stod(values[3]);
				}
				catch (const std::invalid_argument & e) 
				{
					ABORT_MSG(e.what());
				}
				catch (const std::out_of_range & e) 
				{
					ABORT_MSG(e.what());
				}
				
			}
			else
			{
				ABORT_MSG("Can't find simulated annealing parameters");
			}
			
			std::cout<<"Starting simulated annealing..."<<std::endl;
			std::cout<<"Parameters:"<<std::endl;
			std::cout<<"Thread count: "<<info.thread_count<<std::endl;
			std::cout<<"Max outer loops in thread: "<<info.max_outer_loops<<std::endl;
			std::cout<<"Max inner loops in thread: "<<info.max_inner_loops<<std::endl;
			std::cout<<"Initial temperature: "<<info.initial_temperature<<std::endl;
			std::cout<<"Alpha parameter: "<<info.alpha_parameter<<std::endl;
			std::cout<<"Max frozen loops: "<<info.max_frozen_outer_loops<<std::endl;
			std::cout<<"Log level: "<<info.is_log_enabled<<std::endl;
			std::cout<<"Cost Function: "<<info.cost_data->name()<<std::endl;

			if(info.properties_config & SBGEN_USE_NONLINEARITY_FLAG)
			{
				std::cout<<"target NL: "
					<<info.target_properties[SBGEN_NONLINEARITY]
					<<std::endl;
			} else std::cout<<"NL not used\n";

			if(info.properties_config & SBGEN_USE_DELTA_UNIFORMITY_FLAG)
			{
				std::cout<<"target DU: "
					<<info.target_properties[SBGEN_DELTA_UNIFORMITY]
					<<std::endl;
			}else std::cout<<"DU not used\n";

			if(info.properties_config & SBGEN_USE_ALGEBRAIC_IMMUNITY_FLAG)
			{
				std::cout<<"target AI: "
					<<info.target_properties[SBGEN_ALGEBRAIC_IMMUNITY]
					<<std::endl;
			}else std::cout<<"AI not used\n";

			if(info.use_random_seed==true)
				std::cout<<"Seed: random"<<std::endl;
			else
				std::cout<<"Seed: "<<info.seed<<std::endl;
			
			auto sbox = sbgen::simulated_annealing<T>(info);

			if (sbox.has_value()) 
			{
				auto sb = sbox.value();
				
				if(erase_points) 
				{
					sbgen::transform_utils::erase_fixed_points(sb, info.seed);
				}
				
				PRINT_SBOX(sb);
				
				std::cout<<"NL= "
					<<sbgen::properties::nonlinearity(sb)<<std::endl;
				std::cout<<"DU= "
					<<sbgen::properties::delta_uniformity(sb)<<std::endl;
				std::cout<<"AI= "
					<<sbgen::properties::algebraic_immunity(sb)<<std::endl;
				std::cout<<"Fixed Points= "
					<<sbgen::properties::fixed_points(sb)<<std::endl;
				
			} else {
				ABORT_MSG("SBox not found. Try another parameters.");
			}
		}
		else if(parameter->second == "genetic")
		{
			auto cost = get_cost_function<T>(options);
			sbgen::genetic_info_t<T> info;
			info.thread_count = thread_count;
			info.iterations_count = try_per_thread;
			info.is_log_enabled = visibility;
			info.default_log_output = info.is_log_enabled;
			info.delete_parents = false;

			info.selection_method =
				get_selection_method(info, options);
			info.cost_function = cost.first;
			info.cost_data.reset(cost.second);
			setup_properties(
				static_cast<sbgen::properties_info_t*>(&info), options);
			
			parameter = options.find(method_params_flag);
			if (parameter != options.end()) 
			{
				std::vector<std::string> values =
					split(parameter->second, ',');
				if(values.size() != 4)
					ABORT_MSG("Invalid simulated annealing parameters.");

				try 
				{
					info.initial_population_count = std::stoi(values[0]);
					info.mutants_per_parent = std::stoi(values[1]);
					info.selection_count = std::stoi(values[2]);
					info.use_crossover = std::stoi(values[3]);
				}
				catch (const std::invalid_argument & e) 
				{
					ABORT_MSG(e.what());
				}
				catch (const std::out_of_range & e) 
				{
					ABORT_MSG(e.what());
				}

			}
			else
			{
				ABORT_MSG("Can't find simulated annealing parameters");
			}
			setup_crossover_properties(info, options);
			
			std::cout<<"Starting genetic method..."<<std::endl;
		
			std::cout<<"Parameters:"<<std::endl;
			std::cout<<"Thread count: "<<info.thread_count<<std::endl;
			std::cout<<"Mutants per parent: "<<info.mutants_per_parent<<std::endl;
			std::cout<<"Selection count: "<<info.selection_count<<std::endl;
			std::cout<<"Iterations count: "<<info.iterations_count<<std::endl;
			std::cout<<"Init s-box count: "<<info.initial_population_count<<std::endl;
			if (info.use_crossover)
			{
				std::cout<<"Child count: "<<info.child_per_parent<<std::endl;
				std::cout<<"Crossover count: "<<info.crossover_count<<std::endl;
			}
			std::cout<<"Log level: "<<info.is_log_enabled<<std::endl;
			std::cout<<"Cost Function: "<<info.cost_data->name()<<std::endl;
			if(info.properties_config & SBGEN_USE_NONLINEARITY_FLAG)
			{
				std::cout<<"target NL: "
					<<info.target_properties[SBGEN_NONLINEARITY]
					<<std::endl;
			} else std::cout<<"NL not used\n";
			if(info.properties_config & SBGEN_USE_DELTA_UNIFORMITY_FLAG)
			{
				std::cout<<"target DU: "
					<<info.target_properties[SBGEN_DELTA_UNIFORMITY]
					<<std::endl;
			}else std::cout<<"DU not used\n";
			if(info.properties_config & SBGEN_USE_ALGEBRAIC_IMMUNITY_FLAG)
			{
				std::cout<<"target AI: "
					<<info.target_properties[SBGEN_ALGEBRAIC_IMMUNITY]
					<<std::endl;
			}else std::cout<<"AI not used\n";
			
			auto sbox = sbgen::genetic<T>(info);

			if (sbox.has_value()) 
			{
				auto sb = sbox.value();
				
				if(erase_points) 
				{
					sbgen::transform_utils::erase_fixed_points(sb, info.seed);
				}
				
				PRINT_SBOX(sb);
				
				std::cout<<"NL= "
					<<sbgen::properties::nonlinearity(sb)<<std::endl;
				std::cout<<"DU= "
					<<sbgen::properties::delta_uniformity(sb)<<std::endl;
				std::cout<<"AI= "
					<<sbgen::properties::algebraic_immunity(sb)<<std::endl;
				std::cout<<"Fixed Points= "
					<<sbgen::properties::fixed_points(sb)<<std::endl;
			} else {
				std::cout<<"SBox not found. Try another parameters.\n";
			}
		}
	}
}

//=============================================================================

int main(
	int				argc,
	char**			argv) 
{
	options_t		options;
	int				sbox_count = default_sbox_count;

	parse_options(argc, argv, options);

	auto sbox_count_parameter = options.find(sbox_count_flag);
	if (sbox_count_parameter != options.end())
	{
		try 
		{
			sbox_count = std::stoi(sbox_count_parameter->second);
		}
		catch (const std::invalid_argument & e) 
		{
			ABORT_MSG(e.what());
		}
		catch (const std::out_of_range & e) 
		{
			ABORT_MSG(e.what());
		}
	}

	auto file_redirect = options.find(to_file_flag);
	if (file_redirect != options.end())
	{
		auto handle = freopen(file_redirect->second.c_str(),"w",stdout);
		if (handle)
			handle = NULL;
	}

	auto parameter = options.find(help_flag);
	if (parameter != options.end() || argc == 1) 
	{
		print_help();
		return 0;
	}
	
	parameter = options.find(version_flag);
	if (parameter != options.end()) 
	{
		std::cout<<"SBGen "<<SBGEN_VERSION<<std::endl;
		return 0;
	}
	
	parameter = options.find(cost_type_flag);
	if (parameter != options.end()) 
	{
		if (parameter->second == "double")
		{
			for (int i = 0; i < sbox_count; i++)
				run_generator<double>(options);
		}
		else if (parameter->second == "int64_t")
		{
			for (int i = 0; i < sbox_count; i++)
				run_generator<int64_t>(options);
		}
		else
			ABORT_MSG("Unknown cost type. Possible values: double,int64_t");
	} else
	{
		for (int i = 0; i < sbox_count; i++)
			run_generator<double>(options);
	}
}

//=============================================================================