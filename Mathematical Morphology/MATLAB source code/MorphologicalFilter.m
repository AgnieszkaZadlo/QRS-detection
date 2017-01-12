function signal_filtered = MorphologicalFilter(signal,fs)
%UNTITLED9 Summary of this function goes here
% signal - sygna³ wejœciowy
% fs - czêstotliwoœæ próbkowania
% signal_filtered - przefiltrowany sygna³

%Baseline correction
Lo = 0.2*fs;
Lc = 1.5*Lo;
Bo(1,1:Lo) = 1;
Bc(1,1:Lc) = 1;

open1 = dilatation(erosion(signal,Bo),Bo);
close1 = erosion(dilatation(open1,Bc),Bc);

Fbc = signal - close1';

%Noise suppresion
% B_pair = {B1,B2}, diffrent shape, same length
% B1 - triangular, B2 = line

B1=[0 1 5 1 0];
B2=[0 0 0 0 0];

open2 = dilatation(erosion(Fbc,B1),B2);
close2 = erosion(dilatation(Fbc,B1),B2);
signal_filtered = 0.5*(open2+close2);


end

