function [idx_neg,idx_pos]=zero_crossing(der,th)
%ZERO_CROSSING finds indices of zero-crossing moments
%
% Inputs
% ------
%       der - one-dimensional vector with signal (eg.derivative)
%       th  - threshold used to zero tolerance estimation
%
% Outputs
% ------
%       idx_neg - indices of zero-crossing with negative slope
%       idx_pos - indices of zero-crossing with positive slope
%
% Dependencies
% ------------
%
%   This function calls:
%
%   This function is called by:
%       - dpi_based_qrs_detector
%
%
% Authors:    P. Wegrzynowicz, A. Zadlo
% Copyright:  .
% Date:       2016-12-02
idx= (1:numel(der)-1)
sign = der(1:end-1).*der(2:end);
moments=sign<=0 ;
positive = der > th;
negative = der < -th;
idx_pos = moments.*positive(1:end-1).*idx;
idx_neg = moments.*negative(1:end-1).*idx;
idx_pos = nonzeros(idx_pos)';
idx_neg = nonzeros(idx_neg)';
end