#include <iostream>
#include <functional>
#include <string>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <map>


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


const struct option longopts[] =
{
	{"version",						no_argument,			0,	version_flag},
	{"visibility",						no_argument,			0,	visibility_flag},
	{"help",							no_argument,			0,	help_flag},
	{"method",						required_argument,	0,	method_flag},
	{"cost_function",				required_argument,	0,	cost_function_flag},
	{"thread_count",				required_argument,	0,	thread_count_flag},
	{"cost_type",					required_argument,	0,	cost_type_flag},
	{"try_per_thread",			required_argument,	0,	try_per_thread_flag},
	{"method_params",			required_argument,	0,	method_params_flag},
	{"cost_function_params",required_argument,	0,	cost_function_params_flag},
	{0,0,0,0},
};

void print_help()
{
	cout << "Usage: sbgen --method [METHOD] [OPTIONS] \n"
		<< "List of options:\n\n"
			<< "  --visibility\n"
			<< "       Enable verbose mode\n"
			<< "  --version\n"
			<< "       Print version info\n"
			<< "  --help\n"
			<< "       Print help message\n"
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
			<< "       max iterations count without any chages\n"
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
			<< "       Example: --method_params=\"{10, 10000, 1000, 0.99}\"\n"
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
			<< "       Has no free options\n";
			<< "       Please refer to https://github.com/KandiyIIT/SBGen/README.md for more information.\n";
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

template<typename T>
std::function<cost_info_t<T>(cost_function_data_t*, std::array<uint8_t, 256>)> 
get_cost_function(std::map<int, std::string>& options)
{
	auto parameter = options.find(method_flag);
	if (parameter == options.end()) 
	{
		if(parameter->second == "whs")
		{
			
		} 
		else if(parameter->second == "pcf")
		{
			
		}
		else if(parameter->second == "wcf")
		{
			
		}
		ABORT_MSG("Unknown cost function. See help for avaliable cost functions");
	}
	return sbgen::whs<T>;
}

template<typename T>
void run_generator(std::map<int, std::string>& options)
{
	auto parameter = options.find(method_flag);
	if (parameter == options.end()) 
	{
		if(parameter->second == "hill_climbing")
		{
			
		} 
		else if(parameter->second == "simulated_annealing")
		{
			
		}
	}
}

int main(int argc, char** argv) 
{
	std::map<int, std::string> options;
	
	parse_options(argc,argv, options);
	
	auto parameter = options.find(help_flag);
	if (parameter != options.end()) 
	{
		print_help();
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
