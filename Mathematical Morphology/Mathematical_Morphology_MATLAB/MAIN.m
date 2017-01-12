fs = 360;

signal_struct = load('sign_100.mat');
val = signal_struct.val;
MLII = (val(1, :)/200)';
V2 = (val(2, :)/200)';
signal = MLII;

time = 0:1/fs:length(signal)*(1/fs)-(1/fs);

second = MorphologicalFilter(signal, fs);

figure,
subplot(2, 1, 1)
plot(time(1:2160), signal(1:2160))
title('Sygna³ rzeczywisty');
xlabel('Czas [s]');
ylabel('Amplituda [mV]');
subplot(2, 1, 2)
plot(time(1:2160), second(1:2160))
title('Sygna³ po filtracji filtrem dolno- i górnoprzepustowym');
xlabel('Czas [s]');
ylabel('Amplituda [mV]');

el_amp = max(second(1:720)) - min(second(1:720));

v1 = interp1([0 15],[0 el_amp], 1:15);
v2 = interp1([0 15],[el_amp 0], 1:15);
element = [v1 v2(1:end-1)];
t=0:(1/fs):length(element)*(1/fs)-(1/fs);
t=t*1000;
figure,
plot(t, element)
title('Element strukturalny u¿ywany do ekstrakcji cech sygna³u')
xlabel('Czas [ms]')
ylabel('Amplituda [mV]')

after_erosion = erosion(second, element);
after_dilatation = dilatation(second, element);

opening = dilatation(after_erosion, element);
closing = erosion(after_dilatation, element);

fun_first = second - ((opening + closing)/2);
fun = fun_first;
zalamekQ = [];
zalamekR = [];
zalamekS = [];
vect = [];
counter = 0;

for i=1:length(fun)
    
    if(fun(i) ~= 0)
        vect=[vect fun(i)];
        counter = 0;
        
    else if(fun(i) == 0 && counter>3 && length(vect)>22)
            [pks,locs] = findpeaks(vect, 'MinPeakDistance', 10, 'MinPeakHeight', 0.00001);
            N = length(vect);
            if(length(locs)==1)
                max_position = i-N+locs-1;
                zalamekR = [zalamekR max_position];
                potentialQ = fun_first(i-N:max_position);
                [value, position] = min(potentialQ);
                min_posQ = i-N+position-1;
                zalamekQ = [zalamekQ min_posQ];
                potentialS = fun_first(max_position:i-counter);
                [value, position] = min(potentialS);
                min_posS = i-length(potentialS)+position-counter;
                zalamekS = [zalamekS min_posS];
            else if(locs>0)
                    zalamekQ = [zalamekQ i-N+locs(1)-1];
                    zalamekS = [zalamekS i-N+locs(2)-1];
                    potentialR = fun_first(zalamekQ(end):zalamekS(end));
                    [value, position] = min(potentialR);
                    zalamekR = [zalamekR zalamekQ(end)+position-1];
                    end
            end
            vect = [];
        else if(counter<=3 && ~isempty(vect))
            counter = counter+1;
            vect=[vect 0];
            else
                vect = [];
            end
        end
    end
end

t=0:1/fs:length(signal)*(1/fs)-(1/fs);
figure,
subplot(2,1,1)
plot(t(1:4000), signal(1:4000))
title('Sygna³ 100\_MLII z bazy MIT-BIH Arrhythmia Database')
xlabel('Czas [s]')
xlim([0 8])
ylabel('Amplituda [mV]')
subplot(2,1,2)
plot(t(1:4000), fun_first(1:4000))
title('Sygna³ 100\_MLII po przeprowadzeniu operacji morfologicznych')
xlabel('Czas [s]')
xlim([0 8])
ylabel('Amplituda [mV]')

figure,
plot(t(1:4000), signal(1:4000))
title('Wyniki detekcji dla sygna³u 100\_MLII (MIT-BIH Arrhythmia Database)')
xlabel('Czas [s]')
xlim([0 8])
ylabel('Amplituda [mV]')
hold on
scatter(t(zalamekQ), signal(zalamekQ))
hold on
scatter(t(zalamekR), signal(zalamekR))
hold on
scatter(t(zalamekS), signal(zalamekS))
legend('sygna³ 100\_MLII', 'za³amek Q', 'za³amek R', 'za³amek S')