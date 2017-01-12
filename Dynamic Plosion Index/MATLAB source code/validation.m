function [acc,sen] = validation(r_ind_db,r,wnd)
%VALIDATION performs simple validation prcess using two sets of
%annotations. Computes accuracy and sensitivity with the selected tolerance
%window
%
% Inputs
% ------
%       r_ind_db - indices of true samples (annotation form database)
%       r - computed indices of qrs complexes
%       wnd - tolerance window in samples
%
% Outputs
% ------
%       acc - accuracy
%       sen - sensitivity
%
% Dependencies
% ------------
%
%   This function calls:
%
%   This function is called by:
%       - analysis
%
%
% Authors:    P. Wegrzynowicz, A. Zadlo
% Copyright:  .
% Date:       2016-12-02

    if size(r,1)<size(r,2)
        r = r';
    end
    if size(r_ind_db,1)<size(r_ind_db,2)
        r_ind_db = r_ind_db';
    end

    tp = 0;
    fp = 0;
    fn = 0;
    a = ones(numel(r),1)*(-wnd:1:wnd);
    r_wide = a + repmat(r,[1,2*wnd+1]);
    a = ones(numel(r_ind_db),1)*(-wnd:1:wnd);
    r_ind_db_wide = a + repmat(r_ind_db,[1,2*wnd+1]);
    tp = numel(intersect(r_ind_db,r_wide(:)));
    fn = numel(r_ind_db) - tp;
    fp = numel(r) - numel(intersect(r,r_ind_db_wide));
    acc = tp/(tp+fp+fn)*100;
    sen = tp/(tp+fn)*100;
    display(['Accuracy: ',num2str(acc),'%   sensitivity: ',...
        num2str(sen),'%; window=',num2str(wnd),' samples'])
end

