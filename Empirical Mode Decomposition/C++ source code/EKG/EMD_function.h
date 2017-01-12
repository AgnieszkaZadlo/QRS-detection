#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <windows.h>
#include <math.h>
#include "spline.h"


void emd_run (std::vector<double>* ptrSignal, int nIMF);

std::vector<double> singleIMF (std::vector<double> CurrentSignal);

std::vector<double> FirstDerv (std::vector<double> CurrentSignal);

std::vector<double> FindMin (std::vector<double> CurrentSignal);

std::vector<double> FindMax (std::vector<double> CurrentSignal);

std::vector<double> Integration (std::vector<double> IMFArray);

std::vector<int> FindPeaks(std::vector<double> CurrentSignal);