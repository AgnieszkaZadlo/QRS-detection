#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <math.h>

using namespace std;
const double pi = 3.141592653589793238462643383279;

int main () {
  string line;
  string plikW;
  int a,b;
  int k=0;
  int j=0;
  int iloscSygnalow=0;
  int pozycja;
  int wybranySygnal;
  string dtText;
  double dt;
  string * jednostka;
  string * sygnaly;
  vector< vector<double> > vec;
  cout << "Podaj nazwe pliku do wczytania wraz z rozszerzeniem:" << '\n';
  cin >> plikW;

  // Wczytywanie sygna³u pliku csv
  ifstream myfile (plikW.c_str());
  if (myfile.is_open())
  {
      getline (myfile,line);
      for(int i=0;i<line.length();i++){
        if(line[i]==','){
            iloscSygnalow++;
        }

      }
      sygnaly = new string [iloscSygnalow];
      jednostka = new string [iloscSygnalow];
      pozycja=-1;
      j=0;
      for(int i=0;i<line.length();i++){
        if(line[i]== ',' || i == line.length()-1){
                if(i<line.length()-1){
                    j++;
                    if(j>1){
                        sygnaly[j-2]=line.substr(pozycja+1,i-pozycja-1);
                    }
                    pozycja=i;
                }else{
                    j++;
                    sygnaly[j-2]=line.substr(pozycja+1,i-pozycja);
                }
        }
      }
      getline (myfile,line);
      pozycja=-1;
      j=0;
      for(int i=0;i<line.length();i++){
        if(line[i]== ',' || i == line.length()-1){
                if(i<line.length()-1){
                    if(j==0){
                        dtText = line.substr(pozycja+1,i-pozycja-1);
                    }
                    j++;
                    if(j>1){
                        jednostka[j-2]=line.substr(pozycja+1,i-pozycja-1);
                    }
                    pozycja=i;
                }else{
                    j++;
                    jednostka[j-2]=line.substr(pozycja+1,i-pozycja-1);
                }
        }
      }

    for (int i = 0; i < iloscSygnalow; i++) {
        vec.push_back(vector<double>());
    }

    while ( getline (myfile,line) )
    {
        k++;
        j=0;
        pozycja=-1;
        for(int i=0;i<line.length();i++){

            if(line[i]== ',' || i == line.length()-1){
                if(i<line.length()-1){
                    j++;
                    if(j>1){
                       vec[j-2].push_back(atof(line.substr(pozycja+1,i-pozycja-1).c_str()));
                    }
                    pozycja=i;
                }else{
                    j++;
                    vec[j-2].push_back(atof(line.substr(pozycja+1,i-pozycja).c_str()));
                }
            }
        }
    }
    myfile.close();
    }else{
       cout << "Nie mozna otworzyc pliku. Upewnij sie czy podano prawidlowa nazwe pliku wraz z rozszerzeniem";
    }
    cout << "Wczytano dane z pliku" << '\n';

    // Odczytanie okresu próbkowania
    for(int i=1;i<dtText.length();i++){
        if(dtText[i] == '0' || dtText[i] == '1' || dtText[i] == '2' || dtText[i] == '3' || dtText[i] == '4' || dtText[i] == '5'
           || dtText[i] == '6' || dtText[i] == '7' || dtText[i] == '8' || dtText[i] == '9' || dtText[i] == '.'){
            pozycja = i;
        }
    }
    dt = atof(dtText.substr(1,pozycja).c_str());
    cout << "dt: " << dt << '\n';


    double * originalSignal = new double [k];
    double * filteredSignal = new double [k];
    double * signalEnvelope = new double [k];
    cout << "Rozmiar sygnalu:" << k << '\n';
    cout <<"Wybierz ktory sygnal ma zostac poddany filtracji wpisujac liczbe przy nazwie sygnalu:"<< '\n';
    for(int i=0;i<iloscSygnalow;i++){
        cout << i+1 << " - " << sygnaly[i] << '\n';
    }
    cin >> wybranySygnal;
    if(wybranySygnal <iloscSygnalow-1 || wybranySygnal > iloscSygnalow){
        cout << "Podano nieprawidlowa liczbe!" << '\n';
        return 0;
    }
    for(int i=0;i<k;i++){
        originalSignal[i]=vec[wybranySygnal-1][i];
    }

    //Projektowanie filtru
    int M = 37;
    float theta = -pi/4;
    float ox = 0.3;
    float oy = 0.2;
    float wa1 = -0.7;
    float wb1 = 0.7;
    float wa2 = 0.7;
    float wb2= -0.7;
    float lambda = 0.1;
    float w1k;
    float w2l;

    double A = ( pow(cos(theta)/ox,2) + pow((sin(theta)/oy),2) );
    double B = -(sin(2*theta)/pow(ox,2)) + (sin(2*theta)/pow(oy,2));
    double C = ( pow(sin(theta)/ox,2) + pow(cos(theta)/oy,2) );

    double G1 [37][37];
    double G2 [37][37];
    double G [37][37];
    double O [37][37];
    double H [37][37];
    double h [37][37];
    double maxG;

    for(int l=0;l<37;l++){
        for(int m=0;m<37;m++){
            w1k = 2*pi*l/M-pi;
            w2l = 2*pi*m/M-pi;
            G1[l][m]=exp( -(A*pow((w1k-wa1),2) + B*(w1k-wa1)*(w2l-wb1)+ C*pow((w2l-wb1),2) ));
            G2[l][m]=exp( -(A*pow((w1k-wa2),2) + B*(w1k-wa2)*(w2l-wb2)+ C*pow((w2l-wb2),2) ));
            O[l][m]=-((M-1)/2)*w1k-((M-1)/2)*w2l;
            if(l==0 && m==0){
                maxG = G1[l][m];
            }
            if(G1[l][m] > maxG){
                maxG = G1[l][m];
            }
            if(G2[l][m] > maxG){
                maxG = G2[l][m];
            }
        }
    }

    for(int l=0;l<37;l++){
        for(int m=0;m<37;m++){
            G[l][m] = (G1[l][m] + G2[l][m]) / maxG;
        }
    }

    for(int l=0;l<37;l++){
        for(int m=0;m<37;m++){
            H[l][m] = (G[l][m] * cos(O[l][m]));
        }
    }

    // IDFT H w celu otrzymania macierzy wspó³czynników filtra h
    for(int n1=0;n1<37;n1++){
        for(int n2=0;n2<37;n2++){
            h[n1][n2] = 0;
            for(int l=0;l<37;l++){
                for(int m=0;m<37;m++){
                    w1k = 2*pi*l/M-pi;
                    w2l = 2*pi*m/M-pi;
                    h[n1][n2] = h[n1][n2] + H[l][m]*cos(n1*w1k + n2*w2l);
                }
            }
            h[n1][n2] = h[n1][n2] * 1/pow(M,2);
        }
    }

    cout << "Obliczono macierz wspolczynnikow filtra" << '\n';

    //Filtracja sygna³u utworzonym filtrem
    for(int i=37;i<k;i++){
        filteredSignal[i-37]=0;
        for(int l=0;l<37;l++){
            for(int m=0;m<37;m++){
                filteredSignal[i-37] = filteredSignal[i-37] + h[l][m]*originalSignal[i-l]*originalSignal[i-m];
            }
        }
    }
    //Zapisanie do pliku sygna³u po filtracji
    ofstream plik2 (plikW.insert(3,"F").c_str());
    if (plik2.is_open()){
        for(int i=0;i<k-37;i++){
        plik2 << filteredSignal[i] <<"\n";
        }
        plik2.close();
    }
    cout << "Sygnal po filtracji zapisany do pliku "<< plikW   << '\n';

    int maximumY;
    double thv;
    int L;
    int maximumL;
    maximumY = 1;
    //Wyznaczanie progu oraz d³ugoœci okna L
    for(int i=0;i<k-37;i++){
        if(filteredSignal[i]>filteredSignal[maximumY])
            maximumY = i;
    }
    thv = filteredSignal[maximumY]*0.1;
    L=0.120/dt;

    //Obliczanie obwiedni
    for(int i=0;i<k-37;i=i+L){
        if(i+L < k-37){
            maximumL=i;
            for(int l=i;l<i+L;l++){
                if(filteredSignal[l]>filteredSignal[maximumL]){
                    maximumL = l;
                }
            }
            for(int l=i;l<i+L;l++){
                signalEnvelope[l]=filteredSignal[maximumL];
            }
        }else{
            maximumL=i;
            for(int l=i;l<k-37;l++){
                if(filteredSignal[l]>filteredSignal[maximumL]){
                    maximumL = l;
                }
            }
            for(int l=i;l<k-37;l++){
                signalEnvelope[l]=filteredSignal[maximumL];
            }
        }
    }

    //Zapis obwiedni do pliku
    ofstream plik4 (plikW.replace(3,1,"E").c_str());
    if (plik4.is_open()){
        for(int i=0;i<k-37;i++){
        plik4 << signalEnvelope[i] <<"\n";
        }
        plik4.close();
    }
    cout << "Obwiednia sygnalu zapisana do pliku " << plikW   << '\n';
    vector< vector<int> > t; //t[0] = t1; t[1] = t2; t[2] =th;
    for (int i = 0; i < 3; i++) {
        t.push_back(vector<int>());
    }

    //Progowanie sygna³u obwiedni - wyznaczanie czasów pocz¹tkowych i koñcowych
    // t1 i t2
    b=0;
    for(int i=0;i<k-37;i++){
        if(signalEnvelope[i]>thv && a != 1){
            a=1;
            t[0].push_back(i);
        }
        if(  (signalEnvelope[i]<thv && a == 1) || (a == 1 && i == k-37) ){
            t[1].push_back(i-1);
            a=0;
            b++;
        }

    }
    //Wyznaczanie czasu po³ówkowego
   for(int i=0;i<t[0].size();i++){
          t[2].push_back((t[0][i]+t[1][i]+1)/2);
    }
    vector<int> maximumS; // maksimum sygna³u po filtracji
    vector<int> maximumO; // maksimum sygna³u oryginalnego

    //Wyznaczenie za³amków R dla sygna³u po filtracji
     for(int i=0;i<t[2].size();i++){
        if(filteredSignal[t[2][i]]>0){
            maximumS.push_back(t[0][i]);
            for(int l=t[0][i];l<=t[1][i];l++){
                if(filteredSignal[l]> filteredSignal[maximumS[i]]){
                    maximumS[i]= l;
                }
            }
        }else{
            maximumS.push_back(t[0][i]);
            for(int l=t[0][i];l<=t[2][i];l++){
                if(filteredSignal[l]> filteredSignal[maximumS[i]]){
                    maximumS[i]= l;
                }
            }
        }
    }

    //Wyznaczenie za³amków R dla oryginalny sygna³
    for(int i=0;i<maximumS.size();i++){
        maximumS[i] = maximumS[i]+18;
    }

    for(int i=0;i<maximumS.size();i++){
        maximumO.push_back(maximumS[i]);
        for (int l=maximumS[i]-20;l<maximumS[i]+20;l++){
            if(originalSignal[l]>originalSignal[maximumO[i]]){
                maximumO[i] = l;
            }
        }
    }

    //Zapisanie wykrytych za³amków R do pliku
    ofstream plik3 ( plikW.replace(3,1,"R").c_str() );
    if (plik3.is_open()){
        plik3 << "'sample'" <<"\n";
        for(int i=0;i<maximumO.size();i++){
        plik3 << maximumO[i] <<'\n';
        }
        plik3.close();
    }
    cout << "Wykryte zalamki R zapisano w pliku "<< plikW  << '\n';
    cout << "Zakonczono przetwarzanie sygnalu i detekcje zalamkow R";

  cin >> a;
  return 0;
}
