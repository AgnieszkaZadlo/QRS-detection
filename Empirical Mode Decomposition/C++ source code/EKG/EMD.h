#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <windows.h>

class EMD {
	
	public:

	void emd_init(std::vector<double> fptrS, int fnIMF);

	void emd_run ();

	std::vector<double>* singleIMF (std::vector<double>* ptrCurrent);
	 
	std::vector<double>* FirstDerv (std::vector<double>* ptrOrgSignal);

	std::vector<double>* FindMax (std::vector<double>* ptrOrgSignal);

	private:
	std::vector<double>* ptrS;
	std::vector<double>* ptrRes;
	int nIMF;
	int lengthSignal;

	
};