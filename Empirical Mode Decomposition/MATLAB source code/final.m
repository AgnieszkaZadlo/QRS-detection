load('100_MLII.dat','-ascii');
X100_MLII=X100_MLII';
A=X100_MLII(1,1:2000);
load('hp_coef.mat');

%filtracja gornoporzepustowa
A_filtr = filtfilt(hp_coef,1, A);

%wyznaczenie szumów
for i=1:2000
    A_noise(i) = A(i) - A_filtr(i);
end

x1=linspace(0,4.5,1636);
                
                figure(1);
                
                subplot(3,1,1);
                plot(x1,A(1:1636));
                title('Rekord 100,odprowadzenie II');
                ylabel('Amplituda [mV]');
                xlabel('Czas [s]')
                xlim([0 4.5])

                subplot(3,1,2);
                plot(x1,A_filtr(1:1636));
                title('Przefiltrowany sygna³');
                ylabel('Amplituda [mV]');
                xlabel('Czas [s]')
                xlim([0 4.5])

                subplot(3,1,3);
                plot(x1,A_noise(1:1636));
                title('P³ywaj¹ca linia izoelektryczna sygna³u');
                ylabel('Amplituda [mV]');
                xlabel('Czas [s]')
                ylim([-1 1])
                xlim([0 4.5])
 
  [IMFss, resi] = emd2(A_filtr,3);
  
    figure(2);
                subplot(4,1,1);
                plot(x1,A_filtr(1:1636));
                title('Przefiltrowany sygna³ 100 MLII');
                ylabel('Amplituda [mV]');
                xlabel('Czas [s]')
                xlim([0 4.5])

                subplot(4,1,2);
                plot(x1,IMFss(1,1:1636));
                title('IMF1');
                ylabel('Amplituda [mV]');
                xlabel('Czas [s]')
                xlim([0 4.5])

                subplot(4,1,3);
                plot(x1,IMFss(2,1:1636));
                title('IMF2');
                ylabel('Amplituda [mV]');
                xlabel('Czas [s]')
                xlim([0 4.5])

                subplot(4,1,4);
                plot(x1,IMFss(3,1:1636));
                title('IMF3');
                ylabel('Amplituda [mV]');
                xlabel('Czas [s]')
                xlim([0 4.5])

delay =12;
fs = 360;

n = length(IMFss);
w=round(0.1*fs);
window=ones(1,w)/w;
delay=floor(w/2)+delay;

for i=1:2
    for j =3 : n
        if( (IMFss(i,j) * IMFss(i, j-1)  > 0) && (( IMFss(i,j) * IMFss(i, j-2)) > 0 ) ) 
            snt(i, j) = abs( IMFss(i, j) * IMFss(i, j-1) * IMFss(i, j-2));
        else 
            snt(i,j) = 0;
        end
    end
    si(i,:) = conv(snt(i,:),window);
end

%suma
ss=sum(si,1);

%Filtracja dolnoprzepustowa
[d,c]=butter(1,2/180,'low');
ssf=filtfilt(d,c,ss);  
ssf=ssf(delay:end-delay+1); 
                    
                        
                                    figure(3);
                                    x2=linspace(0,4.5,1636);         

                                    subplot(2,1,1);
                                    plot(x2,ss(35:1670));
                                    title('Sygna³ 100 MLII po transformacji nieliniowej i integracji');
                                    ylabel('Amplituda [mV]');
                                    xlabel('Czas [s]')
                                    ylim([0 0.006])
                                    xlim([0 4.5])

                                    subplot(2,1,2)
                                    plot(x2,ssf(35:1670));
                                    title('Sygna³ 100 MLII po filtracji dolnoprzepustowej');
                                    ylabel('Amplituda [mV]');
                                    xlabel('Czas [s]')
                                    ylim([0 0.004])
                                    xlim([0 4.5])

%detekcja szczytów
[Ramp, Rind]=findpeaks(ssf);
                
                
                figure(4)
                plot(A_filtr);
                hold on
                plot(Rind,Ramp,'rs')
                title('Wyniki detekcji pików R dla sygna³u 100 MLII');
                ylabel('Amplituda');
                xlabel('Próbki')

