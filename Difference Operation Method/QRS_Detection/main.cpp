#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <stdlib.h>
#include <math.h>
#include "eigen-eigen-f562a193118d\Eigen\Dense"

using namespace std;
using Eigen::VectorXi;

int main()
{
    ifstream EKGsignal("Signals/100m.csv");
    string line;

    vector <double> Amplitude;
    vector <int> SampleNr;

//save EKGsignal to vectors

    while(getline(EKGsignal,line))
    {
        stringstream  lineStream(line);
        string        cell;
        int column=1;
        while(getline(lineStream,cell,','))
        {
            if(column==1){
                int SampleNrEmement;
                SampleNrEmement = atof(cell.c_str());
                SampleNr.push_back(SampleNrEmement);
            }
            if(column==2){
                double AmplitudeElement;
                AmplitudeElement = atof(cell.c_str());
                Amplitude.push_back(AmplitudeElement);
            }
            column++;
        }
    }

//remove two first elements containing data description

    Amplitude.erase(Amplitude.begin(), Amplitude.begin() + 2);
    SampleNr.erase(SampleNr.begin(), SampleNr.begin() + 2);

//signal differentiation
    int ApmtitudeVectorSize=Amplitude.size();

    vector <double> SignalDif;
    SignalDif.push_back(0);
    for(int ElementNr=1;ElementNr<ApmtitudeVectorSize;ElementNr++){
        double SignalDifElement=Amplitude[ElementNr]-Amplitude[ElementNr-1];
        SignalDif.push_back(SignalDifElement);
    }

//signal thresholding
vector <double> SignalThr;

int SignalDifLength=SignalDif.size();
double PositiveCounter=0;
double NegativeCounter=0;
double PositiveSum=0;
double NegativeSum=0;

for(int SignalDifElementNr=0; SignalDifElementNr<SignalDifLength; SignalDifElementNr++){
    if(SignalDif[SignalDifElementNr]>0){
        PositiveCounter++;
        PositiveSum=PositiveSum+SignalDif[SignalDifElementNr];
    }
    if(SignalDif[SignalDifElementNr]<0){
        NegativeCounter++;
        NegativeSum=NegativeSum+SignalDif[SignalDifElementNr];
    }
}

double MVP=PositiveSum/PositiveCounter;
double MVN=NegativeSum/NegativeCounter;

    double T1=2*MVP;
    double T2=2*MVN;

    for(int AmplitudeVectorElementNr=0; AmplitudeVectorElementNr<ApmtitudeVectorSize; AmplitudeVectorElementNr++){
        if(SignalDif[AmplitudeVectorElementNr]>=T1){
            SignalThr.push_back(SignalDif[AmplitudeVectorElementNr]);
        }
        else if(SignalDif[AmplitudeVectorElementNr]<=T2){
            SignalThr.push_back(SignalDif[AmplitudeVectorElementNr]);
        }
        else{
            SignalThr.push_back(0);
        }
    }

//grouping signal into positive and negative

vector <double> SignalPositive;
vector <double> SignalNegative;

    for(int AmplitudeElementNr=0; AmplitudeElementNr<ApmtitudeVectorSize; AmplitudeElementNr++){
        if(SignalThr[AmplitudeElementNr]>0){
            SignalPositive.push_back(SignalThr[AmplitudeElementNr]);
            SignalNegative.push_back(0);
        }
        else if(SignalThr[AmplitudeElementNr]<0){
            SignalNegative.push_back(SignalThr[AmplitudeElementNr]);
            SignalPositive.push_back(0);
        }
        else{
            SignalNegative.push_back(0);
            SignalPositive.push_back(0);
        }
    }
//maximum and minimum detection in 50 samples intervals
    int AmountOfIntervals=ApmtitudeVectorSize/50;

    vector <double> Maximum;
    vector <double> Minimum;

    vector <double> MaximumNr;
    vector <double> MinimumNr;

    for (int IntervalNr=0; IntervalNr<AmountOfIntervals; IntervalNr++){

        int TempMaximumNr=IntervalNr*50;
        double TempMaximum=SignalPositive[TempMaximumNr];

        for(int IntervalElementNr=0; IntervalElementNr<50; IntervalElementNr++){

            int SignalElement=IntervalNr*50+IntervalElementNr;

            if(SignalPositive[SignalElement]>TempMaximum){
                TempMaximum=SignalPositive[SignalElement];
                TempMaximumNr=SignalElement;
            }
        }

        if(TempMaximum!=0){
            Maximum.push_back(TempMaximum);
            MaximumNr.push_back(TempMaximumNr);
      }

    }

    for (int IntervalNr=0; IntervalNr<AmountOfIntervals; IntervalNr++){

        int TempMinimumNr=IntervalNr*50;
        double TempMinimum=SignalNegative[TempMinimumNr];

        for(int IntervalElementNr=0;IntervalElementNr<50; IntervalElementNr++){

            int SignalElement=IntervalNr*50+IntervalElementNr;

            if(SignalNegative[SignalElement]<TempMinimum){
                TempMinimum=SignalNegative[SignalElement];
                TempMinimumNr=SignalElement;
            }
        }
        if(TempMinimum!=0){
            Minimum.push_back(TempMinimum);
            MinimumNr.push_back(TempMinimumNr);
        }
    }
//check distance between extremes

    vector <double> MaximumTrue;
    vector <double> MinimumTrue;

    vector <double> MaximumTrueNr;
    vector <double> MinimumTrueNr;

    int AmountOfMaximum=Maximum.size();
    int AmountOfMinumum=Minimum.size();

//maximum

    for(int MaximumElementNr=1; MaximumElementNr<AmountOfMaximum; MaximumElementNr++){

        if(MaximumTrueNr.size()==0){
            int MaximumDistance=MaximumNr[MaximumElementNr]-MaximumNr[MaximumElementNr-1];
            if(MaximumDistance>50){
                MaximumTrueNr.push_back(MaximumNr[MaximumElementNr-1]);
                MaximumTrueNr.push_back(MaximumNr[MaximumElementNr]);

                MaximumTrue.push_back(Maximum[MaximumElementNr-1]);
                MaximumTrue.push_back(Maximum[MaximumElementNr]);
            }
            else if(Maximum[MaximumElementNr]>Maximum[MaximumElementNr-1]){
                MaximumTrueNr.push_back(MaximumNr[MaximumElementNr]);
                MaximumTrue.push_back(Maximum[MaximumElementNr]);
            }
            else{
                MaximumTrueNr.push_back(MaximumNr[MaximumElementNr-1]);
                MaximumTrue.push_back(Maximum[MaximumElementNr-1]);
            }
        }
        else{
            int MaximumTrueNrSize=MaximumTrueNr.size();

            int MaximumTrueNrLastElementNr=MaximumTrueNrSize-1;
            double LastMaximum=MaximumTrue[MaximumTrueNrLastElementNr];
            int LastMaximumNr=MaximumTrueNr[MaximumTrueNrLastElementNr];

            int MaximumDistance=MaximumNr[MaximumElementNr]-LastMaximumNr;

            if(MaximumDistance>50){
                MaximumTrueNr.push_back(MaximumNr[MaximumElementNr]);
                MaximumTrue.push_back(Maximum[MaximumElementNr]);
            }
            else if(LastMaximum<Maximum[MaximumElementNr]){
                MaximumTrueNr[MaximumTrueNrLastElementNr]=MaximumNr[MaximumElementNr];
                MaximumTrue[MaximumTrueNrLastElementNr]=Maximum[MaximumElementNr];
            }
        }
    }

//minimum

    for(int MinimumElementNr=1; MinimumElementNr<AmountOfMinumum; MinimumElementNr++){
        if(MinimumTrueNr.size()==0){
            int MinimumDistance=MinimumNr[MinimumElementNr]-MinimumNr[MinimumElementNr-1];

            if(MinimumDistance>50){
                MinimumTrueNr.push_back(MinimumNr[MinimumElementNr-1]);
                MinimumTrueNr.push_back(MinimumNr[MinimumElementNr]);
                MinimumTrue.push_back(Minimum[MinimumElementNr-1]);
                MinimumTrue.push_back(Minimum[MinimumElementNr]);
            }
            else if(Minimum[MinimumElementNr]<Minimum[MinimumElementNr-1]){
                MinimumTrueNr.push_back(MinimumNr[MinimumElementNr]);
                MinimumTrue.push_back(Minimum[MinimumElementNr]);
            }
            else{
                MinimumTrueNr.push_back(MinimumNr[MinimumElementNr-1]);
                MinimumTrue.push_back(Minimum[MinimumElementNr-1]);
            }
        }
        else{
            int MinimumTrueNrSize=MinimumTrueNr.size()-1;
            double LastMinimum=MinimumTrue[MinimumTrueNrSize];
            int LastMinimumNr=MinimumTrueNr[MinimumTrueNrSize];

            int MinimumDistance=MinimumNr[MinimumElementNr]-LastMinimumNr;

            if(MinimumDistance>50){
                MinimumTrueNr.push_back(MinimumNr[MinimumElementNr]);
                MinimumTrue.push_back(Minimum[MinimumElementNr]);
            }
            else if(LastMinimum>Minimum[MinimumElementNr]){
                MinimumTrueNr[MinimumTrueNrSize]=MinimumNr[MinimumElementNr];
                MinimumTrue[MinimumTrueNrSize]=Minimum[MinimumElementNr];

            }
        }
    }

//thresholding of maximum

vector <double> MaximumTresh;
vector <double> MaximumTreshNr;

int AmountOfTrueMaximum=MaximumTrue.size();

double threshold=(*max_element(MaximumTrue.begin(), MaximumTrue.end()))*0.5;

for(int TrueMaximumNr=0; TrueMaximumNr<AmountOfTrueMaximum; TrueMaximumNr++){

    if(MaximumTrue[TrueMaximumNr]>=threshold){
        MaximumTresh.push_back(MaximumTrue[TrueMaximumNr]);
        MaximumTreshNr.push_back(MaximumTrueNr[TrueMaximumNr]);
    }
}

vector <double> detectedComplexes = MaximumTreshNr;

ofstream output_file("Rpoints/100m.atr");
ostream_iterator <int> output_iterator(output_file, "\n");
copy(MaximumTreshNr.begin(), MaximumTreshNr.end(), output_iterator);

// Q and S detection
vector <int> QNumbers;
vector <int> SNumbers;

int RNumber =  MaximumTreshNr.size();

int QStart=0;
int QEnd=(Amplitude.size())-1;

for(int RCounter=1; RCounter<RNumber; RCounter++){
    if(MaximumTreshNr[RCounter]-20>0){
        QStart = MaximumTreshNr[RCounter]-20;
        QEnd = MaximumTreshNr[RCounter];
    }else{
        QStart = 0;
        QEnd = MaximumTreshNr[RCounter];
    }
    double DetectedMin=Amplitude[QStart];

    int DetectedMinNr=QStart;
    for(int QCounter=QStart; QCounter<=QEnd; QCounter++){
        if(Amplitude[QCounter]<DetectedMin){
            DetectedMin=Amplitude[QCounter];
            DetectedMinNr=QCounter;
        }
    }
    QNumbers.push_back(DetectedMinNr);
}
int SStart=0;
int SEnd=(Amplitude.size())-1;
for(int RCounter=0; RCounter<RNumber; RCounter++){
    if(detectedComplexes[RCounter]+20<Amplitude.size()){
        SStart = detectedComplexes[RCounter];
        SEnd = detectedComplexes[RCounter]+20;
    }else{
        SStart = detectedComplexes[RCounter];
        SEnd = (Amplitude.size())-1;
    }
    double DetectedMin=Amplitude[SStart];
    int DetectedMinNr=SStart;
    for(int SCounter=SStart; SCounter<=SEnd; SCounter++){
        if(Amplitude[SCounter]<DetectedMin){
            DetectedMin=Amplitude[SCounter];
            DetectedMinNr=SCounter;
        }
    }
    SNumbers.push_back(DetectedMinNr);
}
}
