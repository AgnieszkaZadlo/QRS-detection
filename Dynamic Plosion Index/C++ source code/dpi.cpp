#include "dpi.h"

VectorXi dpi_based_qrs_detector(VectorXf signal,float fs, float wnd, float p){

    int n285 = floor(0.285*fs);      //The period of the highest possible heart rate (210 BPM)
    int nWnd = floor(wnd/1000.0*fs); //Computation window
    int nWnd285 = n285 + nWnd;       //Computation window with additional 285 ms
    int currQrs = 1,prevQrs = 0;
    int indPrevQrs = 4;
    vector<int> vqrs;
    VectorXf hhecg(nWnd),dpi(nWnd),der(nWnd),filt(nWnd);
//    VectorXi indPos,indNeg;
    MatrixXf dpi_denom(nWnd,nWnd);

    // Highpass filtering in the frequency domain (fc = 8 Hz)
    signal = hpf(signal,FC,fs);

    //Prepare triangle matrix of DPI denominators:
    dpi_denom = getDenominators(nWnd,p);

    while (indPrevQrs < signal.size()-nWnd285){

        VectorXi indPos,indNeg;
        // Half wave of filtered signal
        hhecg = getHalfWaveOfSignal(signal.segment(indPrevQrs-4,nWnd));
        // Dynamic Plosion Index in window
        dpi = (dpi_denom * hhecg).cwiseInverse();
        // Smoothing DPI
        dpi = smooth(dpi);
        // Derivative computation
        der = derivative(dpi);
        tie(indPos,indNeg) = zeroCrossing(der.tail(nWnd-n285), 0.0);
        indPos = indPos.array() + n285;
        indNeg = indNeg.array() + n285;

        int shift = 0;
        if (indPos(0) < indNeg(0)){
            shift++;
        };
        shift += swing(dpi,indPos.segment(shift,indPos.size()-shift),indNeg);
        indPrevQrs += indPos(shift);
        indPrevQrs = improveComplex(indPrevQrs-n285, 2*n285, signal);

        vqrs.push_back(indPrevQrs);
        prevQrs++;
        currQrs++;
    }

    VectorXi qrs(vqrs.size());
    qrs = VectorXi::Map(vqrs.data(),vqrs.size());
    return qrs;
}

tuple<VectorXf,VectorXf,VectorXf> readRecording(string pathToEcgFile){
    string line;
    float col1,col2,col3;
    vector<float> vtime,vlead1,vlead2;
    std::ifstream myfile(pathToEcgFile);
    cout << pathToEcgFile << " : ";
    if (myfile.is_open()){
        getline (myfile,line);
        getline (myfile,line);
        while (getline (myfile,line)){
            std::istringstream ss(line);
            ss >> col1 >> col2 >> col3;
            vtime.push_back(col1);
            vlead1.push_back(col2);
            vlead2.push_back(col3);
            }
        myfile.close();
    }
    else cout << "Unable to open file" << endl;
    VectorXf time(vtime.size()),lead1(vlead1.size()),lead2(vlead1.size());
    time = VectorXf::Map(vtime.data(),vtime.size());
    lead1 = VectorXf::Map(vlead1.data(),vlead1.size());
    lead2 = VectorXf::Map(vlead2.data(),vlead2.size());
    return make_tuple(time,lead1,lead2);
};

VectorXi readAnnotation(string pathToAnnotationFile){
    string line, col1;
    vector<int> ann;
    int col2;
    std::ifstream myfile(pathToAnnotationFile);
    if (myfile.is_open()){
        while (getline (myfile,line)){
            std::istringstream ss(line);
            ss >> col1 >> col2;
            ann.push_back(col2);
        }
        myfile.close();
    }
    VectorXi annotation(ann.size());
    annotation = VectorXi::Map(ann.data(),ann.size());
    return annotation;
}

VectorXf hpf(VectorXf signal, float fc, float fs){

    FFT<float> fft;
    size_t nRows = signal.size();
    float freqPoint = fs/nRows;
    int nPointsToFc = floor(fc/freqPoint);

    VectorXf raisedCosine(nPointsToFc);
    raisedCosine.setLinSpaced(nPointsToFc,1/float(nPointsToFc),1);
    raisedCosine = 0.5-0.5*(raisedCosine.array()*PI).cos();
    VectorXf freqDomainFilter(nRows);
    freqDomainFilter.setOnes(nRows);
    freqDomainFilter.head(nPointsToFc) = raisedCosine;
    freqDomainFilter.tail(nPointsToFc) = raisedCosine.reverse();

    VectorXcf spectrum(nRows);
    fft.fwd(spectrum, signal);
    spectrum = spectrum.cwiseProduct(freqDomainFilter);
    fft.inv(signal,spectrum);
    return signal;
}

VectorXf getHalfWaveOfSignal(VectorXf signal){
    Array<bool, Dynamic,1> result = signal.array()>0.0;
    signal = signal.array().abs() * result.cast<float>();
    signal = signal.array() + 0.0000000001;
    return signal;
}

MatrixXf getDenominators(int wnd, float p){
    VectorXf vec(wnd);
    vec.setLinSpaced(wnd,1.0,float(wnd));
    vec = 1.0/vec.array();
    vec = vec.array().pow(1.0/p);
    MatrixXf dpi_denom(wnd,wnd);
    dpi_denom = vec * RowVectorXf::Ones(wnd);
    dpi_denom = dpi_denom.triangularView<Lower>();
    return dpi_denom;
}

VectorXf convolve(VectorXf vec1, VectorXf vec2, string mode){
    size_t n1 = vec1.size(), n2 = vec2.size();
    VectorXf vec1Wide = VectorXf::Zero(n1+2*n2-2);
    VectorXf conv = VectorXf::Zero(n1+n2-1);
    vec1Wide.head(n2-1).fill(vec1(0));
    vec1Wide.tail(n2-1).fill(vec1(n1-1));
    vec1Wide.segment(n2-1,n1) = vec1;
    for (size_t i = 0; i < (n1+n2-1); i++){
        conv(i) = (vec1Wide.segment(i,n2).cwiseProduct(vec2.reverse())).sum();
    };
    if (mode.compare("same")==0){
        return conv.segment(floor(n2/2),n1);
    }
    else if (mode.compare("valid")==0){
        return conv.segment(n2-1,n1-n2+1);
    }
    else{  // mode: 'full'
        return conv;
    }
}

VectorXf smooth(VectorXf signal){
    VectorXf s(5);
    s << 0.2,0.2,0.2,0.2,0.2;
    return convolve(signal, s, "same");
}

VectorXf derivative(VectorXf signal){
    VectorXf d(5);
    d << -1,0,0,0,1;
    return convolve(signal,d,"same");
}

VectorXi findIndices(VectorXi P){
    VectorXi I;
    I = VectorXi::LinSpaced(P.size(),0,P.size()-1);
    VectorXi IP = I;
    IP.conservativeResize(stable_partition(
      IP.data(),
      IP.data()+IP.size(),
      [&P](int i){return P(i)>0;})-IP.data());
    return IP;
    // Source: http://libigl.github.io/libigl/matlab-to-eigen.html
}

tuple<VectorXi,VectorXi> zeroCrossing(VectorXf der, float threshold){
    size_t n = der.size();
    VectorXi indPos, indNeg, moments(n-1), crossThreshold(n-1);
    VectorXf sign(n-1);
    sign = der.head(n-1).cwiseProduct(der.tail(n-1));
    moments = (sign.array() <= 0).cast<int>();
    crossThreshold = (der.head(n-1).array() > threshold).cast<int>();
    indPos = findIndices(moments.cwiseProduct(crossThreshold));
    crossThreshold = (der.head(n-1).array() < -threshold).cast<int>();
    indNeg = findIndices(moments.cwiseProduct(crossThreshold));
    return make_tuple(indPos,indNeg);
}

int swing(VectorXf dpi, VectorXi indPos, VectorXi indNeg){
    size_t n = min(indPos.size(),indNeg.size());
    VectorXf sw(n);
    sw = (index(dpi,indNeg.head(n)).array() - index(dpi,indPos.head(n)).array()).abs();
    int maxIndex=0;
    sw.maxCoeff(&maxIndex);
    if (indPos(maxIndex)<150){
        int prevIndex = maxIndex+1;
        sw.segment(prevIndex,sw.size()-prevIndex).maxCoeff(&maxIndex);
        maxIndex += prevIndex;
    }
    return maxIndex;
}

VectorXf index(VectorXf vec, VectorXi ind){
    VectorXf subvector(ind.size());
    for (int i=0; i<ind.size();i++){
        subvector(i) = vec(ind(i));
    }
    return subvector;
}

int improveComplex(int indexStart, int nPoints, VectorXf signal){
    int maxIndex=0;
    signal.segment(indexStart,nPoints).cwiseAbs().maxCoeff(&maxIndex);
    return maxIndex+indexStart;
}

tuple<float,float,float> validateDetector(VectorXi trueAnnotation, VectorXi detectedComplexes,int window){
    float sensitivity=0.0, precision=0.0, accuracy=0.0;
    float tp=0,fp=0,fn=0,tn=0;
    for(int i=0; i<trueAnnotation.size(); i++){
        if ((detectedComplexes.array()>=(trueAnnotation(i)-window)).cwiseProduct(detectedComplexes.array()<=(trueAnnotation(i)+window)).any()){
            tp = tp+1;
        }
        else{
            fn = fn+1;
        }
    };
    fp = max(detectedComplexes.size() - tp,float(0));
    accuracy = (tp+tn)/(tp+tn+fp+fn)*100;
    sensitivity = tp/(tp+fn)*100;
    precision = tp/(tp+fp)*100;
    printf("Accuracy: %.2f%% Sensitivity: %.2f%% Precision: %.2f%%  wnd:%d probes",accuracy,sensitivity,precision,window);
    return make_tuple(accuracy,sensitivity,precision);
}
void writeToFile(VectorXi qrs, string pathToResultFile){
    std::ofstream outputFile;
    outputFile.open(pathToResultFile,fstream::out);
    if (outputFile.is_open()){
        outputFile << "Results of QRS detector using Dynamic Plosion Index."<<endl;
        for (int it=0; it<qrs.size(); it++) {
            outputFile << qrs[it]  <<endl;
        }
    outputFile.close();
    }
    else {
         cout << "Unable to open file";
    }
}

