#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <assert.h>
#include <math.h>
#include "eigen-eigen-f562a193118d\Eigen\Dense"

using namespace std;
using Eigen::VectorXi;

vector<float> erosion (vector<float> vect, vector<float> element) {
    int K = vect.size();
    int M = element.size();
    int P = M - 1;

    vector<float> signal_prim(P, vect.back());
    signal_prim.insert(signal_prim.begin(), vect.begin(), vect.end());
    vector<float> signal_out(K, 0);

    for(int k = 0; k<K; k++) {
        float temp[M];
        for(int m = 0; m<M; m++) {
            temp[m] = signal_prim[k+m] - element[m];
        }
        signal_out[k] = *min_element(temp, temp+M);
    }
    return signal_out;
}

vector<float> dilatation (vector<float> vect, vector<float> element) {
    int K = vect.size();
    int M = element.size();
    int P = M - 1;

    vector<float> signal_prim(P, *vect.begin());
    signal_prim.insert(signal_prim.end(), vect.begin(), vect.end());
    vector<float> signal_out(K, 0);

    for(int k = P; k<K+M; k++) {
        float temp[M];
        for(int m = 0; m<M; m++) {
            temp[m] = signal_prim[k-m] + element[m];
        }
        signal_out[k-P] = *max_element(temp, temp+M);
    }
    return signal_out;
}

vector<float> filtration (vector<float> vect){

    //FILTRATION
    //1.Baseline correction

    int Fs = 360;
    int Lo = 0.2*Fs;
    int Lc = 1.5*Lo;

    vector<float> Bo(Lo);
    vector<float> Bc(Lc);

    for(int i=0; i<Lo; i++){
        Bo[i]=1;
    }
    for (int m=0; m<Lc; m++){
        Bc[m]=1;
    }

    vector<float> signal_ero_1 = erosion(vect,Bo);
    vector<float> signal_open_1 = dilatation(signal_ero_1,Bo);
    vector<float> signal_dil_1 = dilatation(signal_open_1,Bc);
    vector<float> signal_close_1 = erosion(signal_dil_1,Bc);
    vector<float> F_base_correct(vect.size());

    for (size_t k=0; k < vect.size(); k++){
        F_base_correct[k] = vect[k] - signal_close_1[k];
    }

    //2.Noise suppression

    vector<float> B1;
    B1.push_back(0);
    B1.push_back(1);
    B1.push_back(5);
    B1.push_back(1);
    B1.push_back(0);

    vector<float> B2;
    B2.push_back(0);
    B2.push_back(0);
    B2.push_back(0);
    B2.push_back(0);
    B2.push_back(0);

    vector<float> signal_ero_2 = erosion(F_base_correct,B1);
    vector<float> signal_open_2 = dilatation(signal_ero_2,B2);
    vector<float> signal_dil_2 = dilatation(F_base_correct,B1);
    vector<float> signal_close_2 = erosion(signal_dil_2,B2);
    vector<float> signal_filtered (vect.size());

    for (size_t j=0; j<vect.size(); j++){
        signal_filtered[j] = 0.5*(signal_open_2[j] + signal_close_2[j]);
    }

    return signal_filtered;

}


vector<float> getElement(float maxValue, int maxPosition) {
    float step = maxValue/maxPosition;
    vector<float> element;
    float value = 0;
    for(int i = 0; i<maxPosition; i++) {
        element.push_back(value);
        value = value + step;
    }
    value = maxValue;
    for(int i = maxPosition; i<2*maxPosition; i++) {
        element.push_back(value);
        value = value - step;
    }
    element.erase(element.begin());
    return element;
}

vector<float> calculateMorphologyFunction(vector<float> signal, vector<float> opening, vector<float> closing) {
    assert(signal.size() == opening.size());
    assert(signal.size() == closing.size());

    int sizeOfSignal = signal.size();
    vector<float> result(sizeOfSignal, 0);
    for(int i = 0; i < sizeOfSignal; i++) {
        result[i] = signal[i] - ((opening[i] + closing[i])/2);
    }
    return result;
}

vector<int> findAllPeaks(vector<float> signal, int minDistance) {
    int vectorSize = signal.size();
    vector<int> output;
    for(int i = 1; i < vectorSize-1; i++) {
        if(signal[i - 1] < signal[i] && signal[i+1] <= signal[i] ) {
            if(output.empty() || (output.back() + minDistance >= i)) {
                output.push_back(i);
            }
        }
    }
    return output;
}

void validateDetector(VectorXi trueAnnotation, VectorXi detectedComplexes,int window){
    float sensitivity=0.0, precision=0.0, accuracy=0.0;
    float tp=0,fp=0,fn=0,tn=0;

    for(int i=0; i<trueAnnotation.size(); i++){
        if ((detectedComplexes.array()>=(trueAnnotation(i)-window)).cwiseProduct(detectedComplexes.array()<=(trueAnnotation(i)+window)).any()){
            tp++;
        }
        else{
            fn++;
        }
    };
    fp = detectedComplexes.size() - tp;
    accuracy = (tp+tn)/(tp+tn+fp+fn)*100;
    sensitivity = tp/(tp+fn)*100;
    precision = tp/(tp+fp)*100;

    printf("Accuracy: %.2f%% Sensitivity: %.2f%% Precision: %.2f%%  wnd:%d probes",accuracy,sensitivity,precision,window);
    return;
}

int main()
{
    fstream plik( "101_MLII.txt", ios::in );
    string val;
    vector <float> vect_1;
    while(!plik.eof()) {
        getline( plik, val );
        if(val != ""){
            vect_1.push_back(atof(val.c_str()));
        }
    }
    plik.close();

    vector <float> vect = filtration(vect_1);

    float element_amplitude = *max_element(vect.begin(), vect.begin() + 720) - *min_element(vect.begin(), vect.begin() + 720);

    vector<float> element = getElement(element_amplitude, 15);

    vector<float> opening = dilatation(erosion(vect, element), element);
    vector<float> closing = erosion(dilatation(vect, element), element);

    vector<float> result = calculateMorphologyFunction(vect, opening, closing);



    vector<int> zalamekQ;
    vector<int> zalamekR;
    vector<int> zalamekS;

    int counter = 0;
    vector<float> temporaryVector;

    for(int i = 0; i<result.size(); i++) {
        if(fabs(result[i])>= 0.0001 ) {
            temporaryVector.push_back(result[i]);
            counter = 0;
        } else if ( fabs(result[i]) < 0.0001 && counter > 3 && temporaryVector.size() > 22) {
            vector<int> peaks = findAllPeaks(temporaryVector, 10);
            int n = temporaryVector.size();
            if(peaks.size() == 1) {
                int maxPosition = i-n+peaks[0];
                zalamekR.push_back(maxPosition);
                //vector<float> potentialQ(vect.begin() + (i - n), maxPosition);
                zalamekQ.push_back(distance(result.begin(),min_element(result.begin() + (i - n), result.begin() + maxPosition)));
                zalamekS.push_back(distance(result.begin(),min_element(result.begin() + maxPosition, result.begin() + (i - counter-1))));
            } else if(peaks.size() > 0) {
                zalamekQ.push_back(i-n+peaks[0]);
                zalamekS.push_back(i-n+peaks[1]);
                zalamekR.push_back(distance(result.begin(), min_element(result.begin() + zalamekQ.back(), result.begin() + zalamekS.back())));
            }
            temporaryVector.clear();
        } else if (counter <= 3 && !temporaryVector.empty()) {
            counter++;
            temporaryVector.push_back(0);
        } else {
            temporaryVector.clear();
        }
    }

    ofstream myfile;
    myfile.open ("zalamekR.txt");
    for(int i = 0; i<zalamekR.size(); i++) {
        myfile<<zalamekR[i]<<endl;
    }
    myfile.close();


    fstream plik1( "atr101.txt", ios::in );
    string val1;
    vector<int> vect1;
    while(!plik1.eof()) {
        getline( plik1, val1 );
        if(val1 != ""){
            vect1.push_back(atoi(val1.c_str()));
        }
    }
    plik1.close();

    Eigen::Map<Eigen::VectorXi> attr(vect1.data(), vect1.size());
    Eigen::Map<Eigen::VectorXi> zalamekRdata(zalamekR.data(), zalamekR.size());

    validateDetector(attr, zalamekRdata, 36);

    system("Pause");

}
