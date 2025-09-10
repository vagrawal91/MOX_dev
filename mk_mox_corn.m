close all;
clear all;

load dataset_corn.mat;
Qlabels = {'Moisture', 'Oil', 'Protein', 'Starch'};

% Vector normalization of NIR spectra
%X = X ./ (std(X')' * ones(1, size(X,2)));

m = size(X, 1);
p = size(X, 2);
q = size(Y, 2);
r = min(p,q);

MSE = nan(r, r);
for k=1:r
  for l=1:min(k, r)
    [P, D, Q, muX, muY, E, Fmax, MSEcv, Fmaxcv, A, B, W] = moxregress(X, Y, k, l, 'CV', 10, 'MCReps', 10);
    MSE(k,l) = MSEcv / q;
  end;
end;


figure(10);
bar3(MSE);
ylabel('{\itk}');
xlabel('{\itl}');
zlabel('MSE (cross-validated)');


k = r;
l = 3;
[P, D, Q, muX, muY, E, Fmax, A, B, W] = mox(X, Y, k, l);


fig1 = figure(1);
[axs, lgd] = moxplot(P, D, Q, X - ones(m, 1)*muX, Y - ones(m, 1)*muY, 'Corn', [], Qlabels, [], [], true, lambda);
xlabel(axs(1), '\lambda [nm]');
xlabel(axs(2), '');

fsize = [12 17]; % cm
set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 fsize], 'PaperSize', fsize);
set(gcf, 'Units', 'centimeters', 'Position', [0 0 fsize]);
%exportgraphics(gcf, 'img/mox_corn.pdf', 'ContentType', 'vector');

