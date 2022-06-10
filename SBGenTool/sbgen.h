#ifndef  _SBGEN_H_
#define _SBGEN_H_

enum method_type_t
{
	METHOD_HILL_CLIMBING,
	METHOD_SIMULATED_ANNEALING
};

#define ABORT_MSG(y)\
  {\
    cerr << "sbgen: " << y << endl;\
    exit(1);\
  }

#endif //_SBGEN_H_
