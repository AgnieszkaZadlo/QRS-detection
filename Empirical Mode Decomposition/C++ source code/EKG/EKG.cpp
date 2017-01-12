// ECG_console.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <windows.h>
#include <iterator>
#include "EMD_function.h"

using namespace std;

std::vector<double> ReadFile(string fileName);

int main(){	
	std::vector<double> bufferData,bufferEMD;
	std::vector<double> *ptrBuffer;
	
	bufferData = ReadFile("SignalCSV/230m.csv");
	
	for(int i=0;i<2000;i++) {
		bufferEMD.push_back(bufferData[i]);
	}

	ptrBuffer = &bufferEMD;
	emd_run(ptrBuffer, 3);
	return 0;
}

std::vector<double> ReadFile(string fileName) {

	string line;
	std::vector<double> buffer;
	ifstream file;
	file.open(fileName);

	if(file.is_open()) {
		while( getline (file, line) ) {
			string value; 
			unsigned first = line.find(',');
			unsigned last = line.find(',', (first+1));
			string strnew = line.substr (first+1,last-first-1);
			double dValue = std::stod(strnew);
			buffer.push_back(dValue);
		}

		buffer.erase(buffer.begin(), buffer.begin()+2);
	}

	std::vector<string> str;
	return buffer;
}

