#include "dpi.h"

int main(int argc,  char** argv){
//    Import signal form text file to MatrixXf,
//    where columns are: time[probes], lead1, lead2
    VectorXf time,lead1,lead2;
    VectorXi ann,qrs;
    string signalPath, annPath, resultPath;

    signalPath = "100.txt";
    annPath = "ann100.txt";
    if (argc>1){
        signalPath = argv[1];
    }
    if (argc>2){
        annPath = argv[2];
    }
    tie(time,lead1,lead2) = readRecording(signalPath);
    qrs = dpi_based_qrs_detector(lead1,FS, 1800.0, 5.0);
//        Get vector of annotations from file
    ann = readAnnotation(annPath);

    float sensitivity,precision,accuracy;
    int wnd = int(FS*0.100);
    tie(accuracy,sensitivity,precision) = validateDetector(ann, qrs, wnd);
    size_t dotPos = annPath.find_last_of(".");
    resultPath = annPath.substr(0, dotPos);
    resultPath += ".dpi";
    writeToFile(qrs , resultPath);
}
