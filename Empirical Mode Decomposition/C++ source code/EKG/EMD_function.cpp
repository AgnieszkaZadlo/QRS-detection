#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <windows.h>
#include <algorithm>
#include <math.h>
#include "EMD_function.h"
#include "persistence1d.hpp"
#include <assert.h>
#include <stdlib.h>

using namespace p1d;

void emd_run (std::vector<double>* ptrSignal, int nIMF) {
	int lengthSignal = ptrSignal->size();
	std::vector<double> ResArray(*ptrSignal);
	int r=0;
	std::vector<std::vector<double>> IMFarray;
	std::vector<double> INTEGRarray0,INTEGRarray1;
	std::vector<int> PeakArray;
	std::vector<double> TempArray;
	
	for(int i=0;i<nIMF;i++) {
		
		long long int il = static_cast<long long int> (i);
		std::string numIMF = std::to_string(il);
		std::string imfNameFile = "IMF" + numIMF + ".mat";
		
		
		ResArray = singleIMF(ResArray);
		IMFarray.push_back(ResArray);

		std::ofstream OutPutFileIMF(imfNameFile);
		std::ostream_iterator<double> output_iterator(OutPutFileIMF, "\n");
		std::copy(ResArray.begin(), ResArray.end(), output_iterator);

	}

	INTEGRarray0 = Integration(IMFarray[1]);

	//delay

	for(r=0;r<36;r++) {
		TempArray.push_back(0);
	}
	for(r=0;r<INTEGRarray0.size();r++) {
		TempArray.push_back(INTEGRarray0[r]);
	}
	
	PeakArray = FindPeaks(TempArray);
	std::ofstream OutPutFilePeakFind("Peaks/230_peak.atr");
	std::ostream_iterator<float> output_iterator_peaks(OutPutFilePeakFind, "\n");
	std::copy(PeakArray.begin(), PeakArray.end(), output_iterator_peaks);
	
	INTEGRarray1 = Integration(IMFarray[1]);

	std::ofstream OutPutFileIntegration("SignalResult/IntegratedSignal.mat");
	if(OutPutFileIntegration.is_open()){
		for(int i=0; i<INTEGRarray0.size(); i++) {
			OutPutFileIntegration << (INTEGRarray0[i] + INTEGRarray1[i]) << std::endl;
		}
	}
	OutPutFileIntegration.close();
}


std::vector<double> singleIMF (std::vector<double> CurrentSignal) {
		
	int max_sift = 15;
	int iter = 0;
	int lengthSignal = CurrentSignal.size();

	//vector do przechowywania kolejnych iteracji imf
	std::vector<double> IMF_iter_array(CurrentSignal);
	std::vector<double> Org_Signal(CurrentSignal);
	std::vector<double> Temp_IMF;

	//vectory do min i max 
	std::vector<double> MinArray;
	std::vector<double> MaxArray;

	//finding extrema
	for(iter=0; iter<max_sift;iter++ ){
		MinArray.clear();
		MaxArray.clear();
		MinArray = FindMin(FirstDerv(IMF_iter_array));
		MaxArray = FindMax(FirstDerv(IMF_iter_array));

		//min array values
		std::vector<double> MinArrayValues;
		for(iter=0; iter<MinArray.size(); iter++) {
			MinArrayValues.push_back(IMF_iter_array[MinArray[iter]]);
		}

		//max array value
		std::vector<double> MaxArrayValues;
		for(iter=0; iter<MaxArray.size(); iter++) {
			MaxArrayValues.push_back(IMF_iter_array[MaxArray[iter]]);
		}

		//cubic interpolation - min envelope
		tk::spline spline_fun_min;
		spline_fun_min.set_points(MinArray, MinArrayValues, true);

		//cubic interpolation - max envelope
		tk::spline spline_fun_max;
		spline_fun_max.set_points(MaxArray, MaxArrayValues, true);

		//vector for average of envelope
		std::ofstream OutPutFileEnv1("SignalResult/Env_low.mat");
		std::ofstream OutPutFileEnv2("SignalResult/Env_high.mat");
		std::vector<double> EnvAverageArray;
		for(iter=0; iter<lengthSignal; iter++) {
			OutPutFileEnv1 << spline_fun_min(iter) << std::endl;
			OutPutFileEnv2 << spline_fun_max(iter) << std::endl;
			EnvAverageArray.push_back((spline_fun_min(iter) + spline_fun_max(iter))/2);
		}

		OutPutFileEnv1.close();
		OutPutFileEnv2.close();

		for(iter=0;iter<lengthSignal;iter++) {
			Temp_IMF.push_back(IMF_iter_array[iter] - EnvAverageArray[iter]);
		}
		
		IMF_iter_array.clear();
		for(iter=0;iter<lengthSignal;iter++) {
			IMF_iter_array.push_back(Temp_IMF[iter]);
		}
	} 
	return IMF_iter_array;
}
	 
std::vector<double> FirstDerv (std::vector<double> CurrentSignal) {

	std::vector<double> FirstDerv;
	int lengthSignal = CurrentSignal.size();
	int i = 0;
	double tmp1, tmp2;

	FirstDerv.push_back(CurrentSignal[1] - CurrentSignal[0]);
		
	for(i=1; i<lengthSignal-2; i++)  {
		tmp1 = CurrentSignal[i]-CurrentSignal[i-1];
		tmp2 = CurrentSignal[i+1]-CurrentSignal[i];
		FirstDerv.push_back((tmp2+tmp1)/2);
	}

	FirstDerv.pop_back();
	FirstDerv.push_back(CurrentSignal[lengthSignal-1] - CurrentSignal[lengthSignal-2]);

	return FirstDerv;
}

std::vector<double> FindMin (std::vector<double> CurrentSignal) {
	std::vector<double> FindMinArray;
	int lengthSignal = CurrentSignal.size();
	int i =0;
	for (i=0;i<lengthSignal-1;i++) {
		if(((CurrentSignal[i]*CurrentSignal[i+1])<0) && (CurrentSignal[i]<0)) {
		FindMinArray.push_back(i);
		}
	}
	return FindMinArray;
}

std::vector<double> FindMax (std::vector<double> CurrentSignal) {
	std::vector<double> FindMaxArray;
	int lengthSignal = CurrentSignal.size();
	int i =0;
	for (i=0;i<lengthSignal-1;i++) {
		if(((CurrentSignal[i]*CurrentSignal[i+1])<0) && (CurrentSignal[i]>0)) {
			FindMaxArray.push_back(i);
		}
 	}
	return FindMaxArray;
}

std::vector<double> Integration (std::vector<double> IMFArray) {
	std::vector<double> SNT;
	std::vector<double> SI;
	int signalLength = IMFArray.size();
	int i,j,k;
	double tmp;
	int imfs = 2; // only two - due to small duration QRS complex
	int fs = 360;
	int w = 0.1 * 360;
	int convLength = signalLength + w - 1;
	std::vector<double> window;
	double delay = floor(static_cast<double> (w/2)) + 12;
	
	int win = 0;
	double g = 0.02777778;
	window.push_back(10);
	window.push_back(0.00034);

	for(win=0; win<w; win++) {
		window.push_back(g);
	}
	
	for(j=imfs; j<signalLength-2; j++) {
		if((((IMFArray[j] * IMFArray[j-1]) > 0) && (IMFArray[j] * IMFArray[j-2]) > 0)) {
			SNT.push_back((IMFArray[j] * IMFArray[j-1] * IMFArray[j-2]));
		} else {
			SNT.push_back(0);
		}
	}

	for (j=(signalLength-2); j<(convLength+2); j++) {
		SNT.push_back(0);
	}
	for(j=36; j<convLength; j++) {
		tmp = 0;
		for(k=0;k<w;k++) {
			tmp = tmp + (SNT[j-k] * window[k]);
		} 
		SI.push_back(tmp);
	}

	std::ofstream OutPutFileSNT("SignalResult/SNT.mat");
	std::ostream_iterator<float> output_iterator_snt(OutPutFileSNT, "\n");
	std::copy(SNT.begin(), SNT.end(), output_iterator_snt);

	std::ofstream OutPutFileSI("SignalResult/SI.mat");
	std::ostream_iterator<float> output_iterator_si(OutPutFileSI, "\n");
	std::copy(SI.begin(), SI.end(), output_iterator_si);

	return SI;
}

std::vector<int> FindPeaks(std::vector<double> CurrentSignal) {
    std::vector<int> FindPeakArray;
	double threshold;
	threshold = *(max_element(CurrentSignal.begin(), CurrentSignal.end()));
	for(int i=1; i < (CurrentSignal.size()-1); i++) {
        if(CurrentSignal[i+1] <= CurrentSignal[i] && CurrentSignal[i-1] < CurrentSignal[i] && CurrentSignal[i] > (threshold/3)) {
                FindPeakArray.push_back(i+1);
        }
    }
    return FindPeakArray;
}