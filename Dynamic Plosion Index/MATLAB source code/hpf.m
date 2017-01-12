function signal_filtrered = hpf(signal,fc,fs)
%HPF function performs hiph-pass filtereing in the freqency domain
%
% Inputs
% ------
%       signal - one-dimensional vector containing signal form one lead
%       fc - cut off frequency (in Hz)
%       fs - sampling frequency (in Hz)
%
% Outputs
% ------
%       signal_filtered - one-dimensional vector containing filtered signal
%
% Dependencies
% ------------
%
%   This function calls:
%       - fft (MATLAB)
%       - ifft (MATLAB)
%
%   This function is called by:
%       - dpi_based_qrs_detector
%
%
% Authors:    P. Wegrzynowicz, A. Zadlo
% Copyright:  .
% Date:       2016-12-02

if mod(numel(signal),2) ~= 0
    signal=signal(1:end-1);end

freq_point= fs/length(signal);
n_points_to_fc= floor(fc/freq_point);

raised_cosine= 0.5-0.5*cos(pi.*(1:n_points_to_fc)/(n_points_to_fc));
freq_domain_filter  = ones(numel(signal),1);
freq_domain_filter(1:n_points_to_fc) = raised_cosine;
freq_domain_filter(end-n_points_to_fc+1:end) = fliplr(raised_cosine);
freq_domain_filter = freq_domain_filter.^4;

spec= fft(signal);
fft_hpf= freq_domain_filter .* spec;
sig_hpf = ifft(fft_hpf);
signal_filtrered= real(sig_hpf);
end