close all;
clear all;

tecator = load('data/tecator.csv');
idx = 1:215;
idx_predictors = 1:100;
idx_responses = 123:125;
X0 = tecator(idx,idx_predictors);
Y0 = tecator(idx,idx_responses);
% Scaling of spectra is not done with zscore. This safegueards against
% wavelengths with small variance
X = X0/std(X0(:));
Y = zscore(Y0);
save('dataset_tecator.mat', 'X', 'Y', '-nocompression');

clear all;
load data/corn.mat;
%X0 = m5spec.data;
%X0 = mp5spec.data;
X0 = mp6spec.data;
Y0 = propvals.data;
% Scaling of spectra is not done with zscore. This safegueards against
% wavelengths with small variance
lambda = 1100:2:2498;
filter_scale = 10.0;
[X_smooth, dX, d2X, lambda] = derivative_spectroscopy(X0, lambda, filter_scale, 3);

X = dX/std(dX(:));
%X = X0/std(X0(:));
Y = zscore(Y0);
save('dataset_corn.mat', 'X', 'Y', 'lambda', '-nocompression');
