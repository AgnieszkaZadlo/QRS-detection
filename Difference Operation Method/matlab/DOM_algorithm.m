close all;
signal = load('208m.mat');
signal = signal.val;
signal = signal./1000;

%notch filtration (fc = 60Hz)
filtr_60=filter(Num60,1,signal);

%highpass filtration (fc = 0.5Hz)
filtr_05=filter(Num05,1,filtr_60);

%morphological filtration
SE = strel('rectangle', [1 3]);
filtr_morph1 = imopen(filtr_05,SE);
filtr_morph = imclose(filtr_morph1,SE);

%signal differentiation
xd = [];
for n=2:1:length(filtr_morph)
    xd(n) = filtr_morph(n)-filtr_morph(n-1);
end

%add vector of zeros at the end of the signal
delay=zeros(1,31);
xd=[xd delay];

%lowpass filtration (fc = 100hz)
xdf=filter(Num100,1,xd);

% after filtration
xdf_new=xdf(32:end);

%Thresholding
positive_counter = 0;
negative_counter = 0;
positive_sum = 0;
negative_sum = 0;
for n=1:1:length(xdf_new)
    if(xdf_new(n)>0)
        positive_counter = positive_counter+1;
        positive_sum = positive_sum + xdf_new(n);
    elseif(xdf_new(n)<0)
        negative_counter = negative_counter+1;
        negative_sum = negative_sum + xdf_new(n);
    end
end

MVP = positive_sum/positive_counter;
MVN = negative_sum/negative_counter;

T1 = MVP*4;
T2 = MVN*4;
 
tresh_signal =[];
for k=1:1:length(xdf_new)
    if(xdf_new(k)>=T1)
        tresh_signal(k)=xdf_new(k);
    elseif(xdf_new(k)<=T2)
        tresh_signal(k)=xdf_new(k);
    else
        tresh_signal(k)=0;
    end
end

%grouping for negative and postive signals
positive_signal = [];
negative_signal = [];

for k=1:1:length(tresh_signal)
    if(tresh_signal(k)>0)
        positive_signal(k)= tresh_signal(k);
        negative_signal(k) = 0;
    elseif(tresh_signal(k)<0)
        positive_signal(k) = 0;
        negative_signal(k) = tresh_signal(k);
    else
        positive_signal(k)=0;
        negative_signal(k)=0;
    end
end

%R detection
maksima = [];
maksima_nr = [];
for i=1:50:length(positive_signal)
    [y,x] = max(positive_signal(i:(i+49)))
    if(y>0)
        maksima = [maksima y];
        maksima_nr = [maksima_nr (x+i)];
    end
end

minima = [];
minima_nr = [];
for i=1:50:length(negative_signal)
    [y,x] = min(negative_signal(i:(i+49)))
    if(y<0)
        minima = [minima y];
        minima_nr = [minima_nr (x+i)];
    end
end

%MAKSIMA
maksima_new = [];
maksima_new_nr = [];
for n=2:1:length(maksima)
    if(length(maksima_new)==0)
        distance = maksima_nr(2)- maksima_nr(1);
        if(distance>50)
            maksima_new_nr = [maksima_nr(1) maksima_nr(2)];
            maksima_new = [maksima(1) maksima(2)];
        elseif(maksima(2)>maksima(1))
            maksima_new_nr = [maksima_new_nr maksima_nr(2)];
            maksima_new = [maksima_new maksima(2)];
        elseif(maksima(1)>maksima(2))
            maksima_new_nr = [maksima_new_nr maksima_nr(1)];
            maksima_new = [maksima_new maksima(1)];
        end
    else
        count=length(maksima_new);
        last_max=maksima_new(count);
        last_max_nr=maksima_new_nr(count);
        
        distance2 = maksima_nr(n)- last_max_nr;
        
        if(distance2>50)
            maksima_new_nr = [maksima_new_nr maksima_nr(n)];
            maksima_new = [maksima_new maksima(n)];
        elseif(maksima(n)>last_max)
            maksima_new_nr(count) = maksima_nr(n);
            maksima_new(count) = maksima(n);    
        end
    end  
end

%MINIMA
minima_new = [];
minima_new_nr = [];
for n=2:1:length(minima)
    if(length(minima_new)==0)
        distance = minima_nr(2)- minima_nr(1);
        if(distance>50)
            minima_new_nr = [minima_nr(1) minima_nr(2)];
            minima_new = [minima(1) minima(2)];
        elseif(minima(2)>minima(1))
            minima_new_nr = [minima_new_nr minima_nr(2)];
            minima_new = [minima_new minima(2)];
        elseif(minima(1)>minima(2))
            minima_new_nr = [minima_new_nr minima_nr(1)];
            minima_new = [minima_new minima(1)];
        end
    else
        count1=length(minima_new);
        last_min=minima_new(count1);
        last_min_nr=minima_new_nr(count1);
        
        distance2 = minima_nr(n)- last_min_nr;
        if(distance2>50)
            minima_new_nr = [minima_new_nr minima_nr(n)];
            minima_new = [minima_new minima(n)];
        elseif(minima(n)<last_min)
            minima_new_nr(count1) = minima_nr(n);
            minima_new(count1) = minima(n);    
        end
    end  
end

R_value = [];
for i=1:1:length(maksima_new_nr)
    R_value(i)=filtr_morph(i);
end

%detekcja Q i S
Q=[];
Q_value = [];
S=[];
S_value = [];
R_number=length(maksima_new_nr);
for i=1:1:R_number
    if(maksima_new_nr(i)-20>0)
        start_Q=maksima_new_nr(i)-20;
        end_Q=maksima_new_nr(i); 
    else
        start_Q=1;
        end_Q=maksima_new_nr(i);
    end
    [qy, qx]=min(filtr_morph(start_Q:end_Q));
    qx=end_Q-qx;
    Q=[Q qx];
    q = filtr_morph(qx);
    Q_value = [Q_value q]; 
    if((length(filtr_morph)-maksima_new_nr(i)>20))
        end_S=maksima_new_nr(i)+20;
        start_S=maksima_new_nr(i);
    else
        end_S=length(filtr_morph);
        start_S=maksima_new_nr(i);
    end
    [sy, sx]=min(filtr_morph(start_S:end_S));
    sx=sx+start_S;
    S=[S sx];
    s = filtr_morph(sx);
    S_value = [S_value s];
end

for n=1:1:length(maksima_new_nr)
    maksima_new_nr(n)=maksima_new_nr(n)+2;
    R_value(n) = filtr_morph(maksima_new_nr(n));
end


