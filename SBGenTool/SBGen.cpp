#include <iostream>
#include <functional>
#include <string>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <map>
#include <sstream> 

#include "hill_climbing.h"
#include "simulated_annealing.h"
#include "cost_function.h"
#include "utils.h"
#include "sbgen_info.h"


#define ABORT_MSG(y)\
	{\
	std::cerr << "sbgen: " << y << std::endl;\
	exit(1);\
	}


#define no_argument				0
#define required_argument		1 
#define optional_argument		2

#define version_flag							0x00
#define help_flag									0x01
#define method_flag							0x02
#define cost_function_flag 					0x03
#define thread_count_flag					0x04
#define cost_type_flag						0x05
#define visibility_flag							0x06
#define try_per_thread_flag				0x07
#define max_frozen_loops_flag			0x08
#define method_params_flag				0x09
#define cost_function_params_flag		0x0a
#define nonlinearity_flag						0x0b
#define delta_uniformity_flag				0x0c
#define algebraic_immunity_flag		0x0d
#define seed_flag								0x0e
#define erase_points_flag					0x0f

#define default_try_per_thread			10000
#define default_thread_count				1
#define default_log								0
#define default_erase_points				0

const struct option longopts[] =
{
	{"version",						no_argument,			0,	version_flag},
	{"visibility",						no_argument,			0,	visibility_flag},
	{"help",							no_argument,			0,	help_flag},
	{"erase_points",				no_argument,	0,	erase_points_flag},
	{"method",						required_argument,	0,	method_flag},
	{"cost_function",				required_argument,	0,	cost_function_flag},
	{"thread_count",				required_argument,	0,	thread_count_flag},
	{"cost_type",					required_argument,	0,	cost_type_flag},
	{"try_per_thread",			required_argument,	0,	try_per_thread_flag},
	{"method_params",			required_argument,	0,	method_params_flag},
	{"cost_function_params",required_argument,	0,	cost_function_params_flag},
	{"nonlinearity",				required_argument,	0,	nonlinearity_flag},
	{"delta_uniformity",			required_argument,	0,	delta_uniformity_flag},
	{"algebraic_immunity",	required_argument,	0,	algebraic_immunity_flag},
	{"seed",							required_argument,	0,	seed_flag},
	{"max_frozen_loops",		required_argument,	0,	max_frozen_loops_flag},

	{0,0,0,0},
};

void print_help()
{
	std::cout << "Usage: sbgen --method [METHOD] [OPTIONS] \n"
		<< "List of options:\n\n"
			<< "  --visibility\n"
			<< "       Enable verbose mode\n"
			<< "  --version\n"
			<< "       Print version info\n"
			<< "  --help\n"
			<< "       Print help message\n"
			<< "  --seed\n"
			<< "       seed for randomness. Warning: in multithread mode there is\n"
			<< "       additional randomnes caused by concurrency\n"
			<< "  --method [hill_climbing|simulated_annealing]\n"
			<< "       hill_climbing = hill climbing method\n"
			<< "       simulated_annealing = simulated annealing method\n"
			<< "  --cost_function [whs|wcf|pcf]\n"
			<< "       whs = WHS cost function\n"
			<< "       wcf = WCF cost function\n"
			<< "       pcf = PCF cost function\n"
			<< "  --thread_count\n"
			<< "       max thread count\n"
			<< "  --cost_type [int64_t|double]\n"
			<< "       type of variable, where stored s-box cost. Default value - double\n"
			<< "  --try_per_thread\n"
			<< "       maximal iterations count in method\n"
			<< "  --max_frozen_loops\n"
			<< "       max iterations count without any chages\n\n"
		<< "Method parameter list:\n\n"
			<< "  --method_params\n"
			<< "       params of method in format --method_params={param1,param2,...,paramN}\n"
			<< "  hill_climbing:\n"
			<< "       Has no free options\n"
			<< "  simulated_annealing\n"
			<< "       param1: max_outer_loops - maximal outer loop count\n"
			<< "       param2: max_inner_loops - maximal inner loop count\n"
			<< "       param3: initial_temperature -  initial temperature\n"
			<< "       param4: alpha_parameter -  alpha_parameter\n"
			<< "       Example: --method_params=\"{10, 10000, 1000, 0.99}\"\n\n"
		<< "Cost function parameter list:\n\n"
			<< "  --cost_function_params\n"
			<< "       params of cost function in format --cost_function_params={param1,param2,...,paramN}\n"
			<< "  whs\n"
			<< "       param1: r\n"
			<< "       param2: x\n"
			<< "       Example: --cost_function_params=\"{12, 0}\"\n"
			<< "  pcf\n"
			<< "       param1: n\n"
			<< "       Example: --cost_function_params=\"{5}\"\n"
			<< "  wcf\n"
			<< "       Has no free options\n\n"
		<< "Target properties:\n\n"
			<< "  --nonlinearity\n"
			<< "       target nonlinearity value.\n"
			<< "  --delta_uniformity\n"
			<< "       target delta uniformity value.\n"
			<< "  --algebraic_immunity\n"
			<< "       target algebraic immunity value.\n\n"
			<< " Please refer to https://github.com/KandiyIIT/SBGen/README.md for more information.\n";
}

void parse_options(int argc, char **argv, std::map<int, std::string>& options)
{
	int flag;
	int option_index;
	const char* short_options = "";
	
	while ((flag=getopt_long(argc,argv,short_options, longopts,&option_index))!=-1)
	{
		std::string parameter_value = "";
		
		if(flag=='?')
			ABORT_MSG("Unknown parameter");
		
		if(longopts[option_index].has_arg == required_argument)
			parameter_value = optarg;
		options[longopts[option_index].val] = parameter_value;
	}
}

std::vector<std::string> &split(const std::string &s, char delim,std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (item.length() > 0) {
            elems.push_back(item);  
        }
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

template<typename T>
std::pair<std::function<sbgen::cost_info_t<T>(sbgen::cost_function_data_t*, std::array<uint8_t, 256>)>,
	sbgen::cost_function_data_t*>
get_cost_function(std::map<int, std::string>& options)
{
	auto parameter = options.find(cost_function_flag);
	if (parameter != options.end()) 
	{
		if(parameter->second == "whs")
		{
			int32_t x = 0;
			int32_t r = 3;
            
			parameter = options.find(cost_function_params_flag);
			if (parameter != options.end()) 
			{
                std::vector<std::string> values = split(parameter->second, ',');
				
				if(values.size()!=2)
				{
					ABORT_MSG("Invalid parameters count for whs cost function");
				}
				
				try {
					r = std::stoi(values[0]);
					x = std::stoi(values[1]);
				}
				catch (const std::invalid_argument & e) 
				{
					std::cout << e.what() << "\n";
				}
				catch (const std::out_of_range & e) 
				{
					std::cout << e.what() << "\n";
				}
				
				return std::make_pair<std::function< sbgen::cost_info_t<T>(sbgen::cost_function_data_t*, std::array<uint8_t, 256>)>
					, sbgen::cost_function_data_t*>(sbgen::whs<T>,new sbgen::whs_function_data_t(r, x));
			} 
			else
			{
				ABORT_MSG("Can't find whs parameters");
			}
		} 
		else if(parameter->second == "pcf")
		{
			int32_t n = 5;
            
			parameter = options.find(cost_function_params_flag);
			if (parameter != options.end()) 
			{
                std::vector<std::string> values = split(parameter->second, ',');
				
				if(values.size()!=1)
				{
					ABORT_MSG("Invalid parameters count for pcf cost function");
				}
				
				try {
					n = std::stoi(values[0]);
				}
				catch (const std::invalid_argument & e) 
				{
					std::cout << e.what() << "\n";
				}
				catch (const std::out_of_range & e) 
				{
					std::cout << e.what() << "\n";
				}
				
				return std::make_pair<std::function< sbgen::cost_info_t<T>(sbgen::cost_function_data_t*, std::array<uint8_t, 256>)>
				, sbgen::cost_function_data_t*>(sbgen::pcf<T>,new sbgen::pcf_function_data_t(n));
			} 
			else
			{
				ABORT_MSG("Can't find pcf parameters");
			}
		}
		else if(parameter->second == "wcf")
		{
			return std::make_pair<std::function< sbgen::cost_info_t<T>(sbgen::cost_function_data_t*, std::array<uint8_t, 256>)>
				, sbgen::cost_function_data_t*>(sbgen::wcf<T>,new sbgen::wcf_function_data_t());
		}
		ABORT_MSG("Unknown cost function. See help for avaliable cost functions");
	}
	return std::make_pair<std::function< sbgen::cost_info_t<T>(sbgen::cost_function_data_t*, std::array<uint8_t, 256>)>
		, sbgen::cost_function_data_t*>(sbgen::wcf<T>,new sbgen::wcf_function_data_t());
}

void setup_properties (sbgen::properties_info_t* info,  std::map<int, std::string>& options) 
{
	info->properties_config = 0;
	try 
	{
		auto property = options.find(nonlinearity_flag);
		if (property != options.end()) 
		{
			setup_property( info, SBGEN_NONLINEARITY, std::stoi(property->second));
		}
		else
		{
			ABORT_MSG("Need target nonlinearity");
		}
		
		property = options.find(delta_uniformity_flag);
		if (property != options.end()) 
		{
			setup_property( info, SBGEN_DELTA_UNIFORMITY, std::stoi(property->second));
		}
		
		property = options.find(algebraic_immunity_flag);
		if (property != options.end()) 
		{
			setup_property( info, SBGEN_ALGEBRAIC_IMMUNITY, std::stoi(property->second));
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
		std::cout << e.what() << "\n";
	}
	catch (const std::out_of_range & e) 
	{
		std::cout << e.what() << "\n";
	}
	return;
}

template<typename T>
void run_generator(std::map<int, std::string>& options)
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
			std::cout << e.what() << "\n";
		}
		catch (const std::out_of_range & e) 
		{
			std::cout << e.what() << "\n";
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
			setup_properties(static_cast<sbgen::properties_info_t*>(&info), options);
			
			std::cout<<"Starting hill climbing..."<<std::endl;
			std::cout<<"Parameters:"<<std::endl;
			std::cout<<"Thread count: "<<info.thread_count<<std::endl;
			std::cout<<"Try per thread: "<<info.try_per_thread<<std::endl;
			std::cout<<"Max frozen loops: "<<info.max_frozen_count<<std::endl;
			std::cout<<"Log level: "<<info.is_log_enabled<<std::endl;
			std::cout<<"Cost Function: "<<info.cost_data->name()<<std::endl;
			if(info.properties_config & SBGEN_USE_NONLINEARITY_FLAG) {
				std::cout<<"target NL: "<<info.target_properties[SBGEN_NONLINEARITY]<<std::endl;
			} else std::cout<<"NL not used\n";
			if(info.properties_config & SBGEN_USE_DELTA_UNIFORMITY_FLAG) {
				std::cout<<"target DU: "<<info.target_properties[SBGEN_DELTA_UNIFORMITY]<<std::endl;
			}else std::cout<<"DU not used\n";
			if(info.properties_config & SBGEN_USE_ALGEBRAIC_IMMUNITY_FLAG) {
				std::cout<<"target AI: "<<info.target_properties[SBGEN_ALGEBRAIC_IMMUNITY]<<std::endl;
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
				
				std::cout<<"NL= "<<sbgen::properties::nonlinearity(sb)<<std::endl;
				std::cout<<"DU= "<<sbgen::properties::delta_uniformity(sb)<<std::endl;
				std::cout<<"AI= "<<sbgen::properties::algebraic_immunity(sb)<<std::endl;
				std::cout<<"Fixed Points= "<<sbgen::properties::fixed_points(sb)<<std::endl;
				
			} else {
				ABORT_MSG("SBox not found. Try another parameters or restart.");
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
			setup_properties(static_cast<sbgen::properties_info_t*>(&info), options);
			
			parameter = options.find(method_params_flag);
			if (parameter != options.end()) 
			{
                std::vector<std::string> values = split(parameter->second, ',');
				if(values.size() != 4)
					ABORT_MSG("Invalid simulated annealing parameters. See help");
				
				try 
				{
					info.max_outer_loops = std::stoi(values[0]);
					info.max_inner_loops = std::stoi(values[1]);
					info.initial_temperature = std::stod(values[2]);
					info.alpha_parameter = std::stod(values[3]);
				}
				catch (const std::invalid_argument & e) 
				{
					std::cout << e.what() << "\n";
				}
				catch (const std::out_of_range & e) 
				{
					std::cout << e.what() << "\n";
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
			if(info.properties_config & SBGEN_USE_NONLINEARITY_FLAG) {
				std::cout<<"target NL: "<<info.target_properties[SBGEN_NONLINEARITY]<<std::endl;
			} else std::cout<<"NL not used\n";
			if(info.properties_config & SBGEN_USE_DELTA_UNIFORMITY_FLAG) {
				std::cout<<"target DU: "<<info.target_properties[SBGEN_DELTA_UNIFORMITY]<<std::endl;
			}else std::cout<<"DU not used\n";
			if(info.properties_config & SBGEN_USE_ALGEBRAIC_IMMUNITY_FLAG) {
				std::cout<<"target AI: "<<info.target_properties[SBGEN_ALGEBRAIC_IMMUNITY]<<std::endl;
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
				
				std::cout<<"NL= "<<sbgen::properties::nonlinearity(sb)<<std::endl;
				std::cout<<"DU= "<<sbgen::properties::delta_uniformity(sb)<<std::endl;
				std::cout<<"AI= "<<sbgen::properties::algebraic_immunity(sb)<<std::endl;
				std::cout<<"Fixed Points= "<<sbgen::properties::fixed_points(sb)<<std::endl;
				
			} else {
				ABORT_MSG("SBox not found. Try another parameters or restart.");
			}
		}
	}
}

int main(int argc, char** argv) 
{
	std::map<int, std::string> options;
	
	parse_options(argc,argv, options);
	
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
			run_generator<double>(options);
		else if (parameter->second == "int64_t")
			run_generator<int64_t>(options);
		else
			ABORT_MSG("Unknown cost type. Possible values: double,int64_t");
	} else
		run_generator<double>(options);
}
