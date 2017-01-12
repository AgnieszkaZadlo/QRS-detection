function [IMF, residual] = emd (signal,numIMF)

L = length(signal);

IMF = zeros(numIMF, L);
residual = signal;

for n =1: numIMF
    [mode, residual] = Single_IMF (residual);
    IMF(n, :) = mode;
end
end

function[mode, residual] = Single_IMF(signal)

L = length(signal);
samples = linspace (1,L,L);
max_sift = 15;
threshold = 0.09;

cross_number =0;

for z = 1: L
    if (signal(z) == 0) 
        cross_number = cross_number +1;
    end
end

IMF_fun = zeros(max_sift+1,L);
IMF_fun(1,:) = signal;

for j = 1:max_sift
    [min_indexes, max_indexes] = extremum(IMF_fun(j,:));    %finding extrema
    min_values = IMF_fun(j,min_indexes);
    max_values = IMF_fun(j,max_indexes);
    env_min = pchip(min_indexes,min_values,samples);    %lower envelope
    env_max = pchip(max_indexes,max_values,samples);    %upper envelope
    env_mean = (env_min + env_max) ./ 2;    %mean envelope
    IMF_fun(j+1,:) = IMF_fun(j,:) - env_mean;
    max_env_mean = abs(max(env_mean));
    min_env_mean = (min(env_mean));
    if and ((max_env_mean < threshold),(min_env_mean< threshold))
        break
    end    
        
end

mode = IMF_fun(end,:);
residual = signal - mode;

end

function [indMin, indMax] = extremum(signal)

d = diff(signal);      
n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
indMin = find(d1.*d2<0 & d1<0)+1;
indMax = find(d1.*d2<0 & d1>0)+1;

imin = [];
imax = [];
    
if any(d==0)
    dZero = (d==0)
    for i=1:n 
        if(d(i) == 0)
            if ((d(i-1)*d(i+1) <0) & d(i-1) < 0)
            imin = [imin i];
            end
            if ((d(i-1)*d(i+1) <0) & d(i-1) > 0)
            imax = [imax i];
            end
        end
    end
end

  if ~isempty(imax)
    indMax = sort([indMax imax]);
  end

  if ~isempty(imin)
    indMin = sort([indMin imin]);
  end
  
end

