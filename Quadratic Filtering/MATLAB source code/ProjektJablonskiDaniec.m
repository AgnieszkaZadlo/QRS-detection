clear all;
load('100_MLII.dat');
signal = X100_MLII(1:2048);

%Projektowanie filtru
M=37;
theta=-pi/4;
ox=0.3;
oy=0.2;

A=(cos(theta)/ox)^2+(sin(theta)/oy)^2;
B=-(sin(2*theta)/ox^2)+sin(2*theta)/oy^2;
C=(sin(theta)/ox)^2+ (cos(theta)/oy)^2;

G1=zeros();
G2=zeros();
O=zeros();
wa1=-0.7;
wb1=0.7;
wa2=0.7;
wb2=-0.7;

for k=0:36
    for l=0:36
        w1k=2*pi*k/M-pi;
        w2l=2*pi*l/M-pi;
        G1(k+1,l+1)=exp(-(A*(w1k-wa1)^2+B*(w1k-wa1)*(w2l-wb1)+C*(w2l-wb1)^2));
        G2(k+1,l+1)=exp(-(A*(w1k-wa2)^2+B*(w1k-wa2)*(w2l-wb2)+C*(w2l-wb2)^2));
        O(k+1,l+1)=-18*w1k-18*w2l;
    end
end
G=G1+G2;
G=G/max(max(G));
H=G.*exp(j*O);

% IDFT H w celu otrzymania macierzy wspó³czynników filtra h
h=zeros();
for n1=0:36
    for n2=0:36
        sum=0;
        for i=0:M-1
            for k=0:M-1
                w1 = ((2*pi*i)/M) - pi;
                w2 = ((2*pi*k)/M) - pi;
                sum = sum + H(i+1,k+1)*exp(n1*w1*j+j*n2*w2);
            end
        end
        sum=sum*(1/M^2);
        h(n1+1,n2+1)=sum;
    end
end

%Filtracja sygna³u
y = zeros();
for i=38:numel(signal)
    y(i-37)=0;
    for j=1:37
        for k=1:37
            y(i-37)=y(i-37)+h(j,k)*signal(i-j)*signal(i-k);
        end
    end
end

%Wyznaczanie progu
lambda=0.1;
ym=max(y);
thv=lambda*ym;
thv=real(thv);

%Obliczanie obwiedni
L=36; 
z=zeros();

for i=L:length(y)
      z(i)=max(y(i-L+1:i));
end

%% 
%Progowanie sygna³u obwiedni - wyznaczanie czasów pocz¹tkowych i koñcowych
j=1;
k=0;
t1=zeros();
t2=zeros();
for i=1:numel(z)
   if(z(i)>thv && k ~= 1)
       k=1;
       t1(j)=i;
   end
   if((z(i)<thv && k == 1) || (k == 1 && i == numel(z)))
       t2(j)=i-1;
       k=0;
       j=j+1;
   end
end

%Wyznaczanie czasów po³ówkowych
th=zeros();
for i=1:numel(t1)
   th(i)=(t1(i)+t2(i))/2; 
end
th = ceil(th);

%%
%Wyznaczanie maksimów na przetworzonym sygnale
maximum=zeros();
for i=1:numel(th)
    if(y(th(i))>0)
        maximum(i)= t1(i);
        for j=t1(i):t2(i)
            if(y(j)> y(maximum(i)))
                maximum(i)= j;
            end
        end
    else
        maximum(i)= t1(i);
        for j=t1(i):th(i)
            if(y(j)> y(maximum(i)))
                maximum(i)= j;
            end
        end      
    end
end

% Wyznaczenie maksimów na oryginalnym sygnale
maxsignal = maximum + 18;
maxostateczne=zeros();
for i=1:numel(maxsignal)
    maxostateczne(i) = maxsignal(i);
    for j = maxsignal(i)-20:maxsignal(i)+20
        if(signal(j)>signal(maxostateczne(i)))
           maxostateczne(i)=j; 
        end
    end
end


figure,plot(signal(1:end-37))
title('Oryginalny sygna³ z wykrytymi za³amkami R');
xlabel('Numer próbki');
ylabel('Napiêcie [mV]');
axis([-inf inf -inf inf])
hold on;
scatter(maxostateczne,signal(maxostateczne))
