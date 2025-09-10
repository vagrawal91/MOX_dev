function MSE_all = MSEvalues_all(d_X, d_XY, MSEmox_mean, MSEmox_kisl_mean, ...
    MSEpls_mean, MSEcca_mean, MSEols_mean)
% computeMSEs - Computes MSE values at d_XY, d_XY + d_X, and local minimum of PLS
%
% Inputs:
%   d_X               - Total independent latent variables in predictor X
%   d_XY              - Total shared latent variables between pred. & resp.
%   MSEmox_mean       - Vector of MSE values for MOX_rh
%   MSEmox_kisl_mean  - Vector of MSE values for MOX_hh
%   MSEpls_mean       - Vector of MSE values for PLS
%   MSEcca_mean       - Vector of MSE values for CCA
%   MSEols_mean       - Vector of MSE values for OLS
%
% Output:
%   MSE_all           - Combined row vector of MSEs at three points with NaN spacers

%--- MSE at d_XY
dd = d_XY;
MSE_dXY = [MSEmox_mean(dd); MSEmox_kisl_mean(dd); MSEpls_mean(dd); MSEcca_mean(dd); MSEols_mean]';
%disp('MSE at d_XY:'); disp(MSE_dXY);

%--- MSE at d_XY + d_X
dd = d_XY + d_X;
MSE_dxXY = [MSEmox_mean(dd); MSEmox_kisl_mean(dd); MSEpls_mean(dd); MSEcca_mean(dd); MSEols_mean]';
%disp('MSE at d_XY + d_X:'); disp(MSE_dxXY);

%--- MSE at local minimum of PLS
dd = find(MSEpls_mean == min(MSEpls_mean));
MSE_lmPLS = [MSEmox_mean(dd); MSEmox_kisl_mean(dd); MSEpls_mean(dd); MSEcca_mean(dd); MSEols_mean]';
%disp('MSE at local minimum of PLS:'); disp(MSE_lmPLS);

%--- Combine all with NaN spacers (NaN to add space in excel)
MSE_all = [MSE_dXY, NaN, MSE_dxXY, NaN, MSE_lmPLS];
end