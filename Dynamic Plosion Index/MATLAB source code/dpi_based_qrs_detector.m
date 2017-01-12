function qrs = dpi_based_qrs_detector(signal,fs,window,p,figures)
%DPI_BASED_QRS_DETECTOR function find QRS complexes in ECG signal and
%returns indices of them
%
% Inputs
% ------
%       signal - one-dimensional vector containing signal form one ECG lead
%       fs - sampling frequency (in Hz)
%       window - computation window fo DPI algorithm
%       p - dynamic scaling value for DPI algorithm
%       figures - decision to plot steps of detection (bool)
%
% Outputs
% ------
%       qrs - indices of QRS complexes in ECG signal
%
% Dependencies
% ------------
%
%   This function calls:
%       - hpf
%       - zero_crossing
%       - smooth (MATLAB)
%       - conv (MATLAB)
%
%   This function is called by:
%
%
% Authors:    P. Wegrzynowicz, A. Zadlo
% Copyright:  .
% Date:       2016-12-02
if size(signal,1) == 1
    signal=signal';
end
if mod(length(signal),2) ~= 0
    signal=signal(1:end-1);
end

%% Highpass filtering in the frequency domain (fc = 8 Hz)
signal_filt = hpf(signal,8,fs);
% Highpass filtering in the frequency domain (fc = 2 Hz)
signal_filt_2Hz = hpf(signal,2,fs);
%% DPI
fs_kHz=fs/1000;                         %Sampling frequency in kHz
ms_285 =floor(285*fs_kHz);              %The period of the highest possible heart rate (210 BPM)
ms_wnd_285=floor((window+285)*fs_kHz);  %Computation window with additional 285 ms
ms_wnd=floor(window*fs_kHz);            %Computation window

%Triangle matrix of DPI denominators:
vec = 1./(1:ms_wnd+1).^(1/p);
dpi_denom = repmat(vec,[ms_wnd+1,1])';
dpi_denom = tril(dpi_denom);


qrs(1) = 5;
m=2;
i=1;
while i< length(signal_filt)-ms_wnd_285
    %% Cut the part of the signal
    sig_part_f8hz= signal_filt(qrs(m-1)-4:qrs(m-1)-4+ms_wnd);
    
    %% Half wave
    hhecg = sig_part_f8hz.*(sig_part_f8hz>=0);
    
    %% Dynamic Plosion Index
    dpi= 1./(dpi_denom * hhecg);
    dpi=smooth(dpi);
    
    %% DPI derivative
    der = conv(dpi,[-1 0 0 0 1],'same')';
    
    %% Zero crossing points
    [idx_neg,idx_pos]= zero_crossing(der(ms_285:end),0);
    % Shift inidces
    idx_pos=idx_pos+ms_285;
    idx_neg=idx_neg+ms_285;
    
    % Remove first zero-crossing index if it is positive
    % (current R peak infuence)
    if idx_pos(1) < idx_neg(1)
        idx_pos(1)=[];
    end
    
    %% Finding pairs
    n_pairs = min(numel(idx_pos),numel(idx_neg));
    % Find pair with maximum swing value
    swing = abs(dpi(idx_neg(1:n_pairs)) - dpi(idx_pos(1:n_pairs)));
    [~,order] = sort(swing,'descend');
    idx_pos_sort = idx_pos(order);
    idx_neg_sort = idx_neg(order);
    idx_pos_sort = idx_pos_sort(idx_pos_sort>ms_285+1);
    if ~isempty(idx_pos_sort)
        idx = idx_pos_sort(1);
    else
        idx=200;
        warning('QRS was artificially shifted - 200 samples')
    end
    
    %% Absolute indices
    idx = idx+qrs(m-1);
    s = max([1,idx-ms_285]);
    e = min([idx+ms_285,numel(signal_filt_2Hz)]);
    [~,shift] = max(abs(signal_filt_2Hz(s:e)));
    qrs(m)= s+shift;
    i = qrs(m);
    m=m+1;
    
    if figures
        sig_part_orig = signal(qrs(m-1)-4:qrs(m-1)-4+ms_wnd);
        figure(1);
        clf()
        subplot(4,1,1)
        plot(sig_part_orig,'k')
        hold on
        plot(sig_part_f8hz,'g')
        plot(hhecg,'r')
        grid on
        axis tight
        
        subplot(4,1,2)
        plot(dpi)
        axis tight
        ylim([dpi(end),max(dpi(30:end))])
        subplot(4,1,3)
        plot(der)
        hold on
        plot(der(1:ms_285),'r','LineWidth',1.5)
        plot(idx_neg,der(idx_neg),'bo')
        plot(idx_pos,der(idx_pos),'ro')
        axis tight
        ylim([min(der(30:end)),max(der(30:end))])
        
        subplot(4,1,4)
        plot(dpi)
        hold on
        plot(dpi(1:ms_285),'r','LineWidth',1.5)
        axis tight
        ylim([dpi(end),max(dpi(30:end))])
        
        plot(idx_neg_sort(1),dpi(idx_neg_sort(1)),'bo')
        plot(idx_pos_sort(1),dpi(idx_pos_sort(1)),'ro')
        pause();
    end
    
end
% Delete first index, because it was artificial:
qrs(1)=[];
end

