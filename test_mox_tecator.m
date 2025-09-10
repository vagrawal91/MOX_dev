close all;
clear all;

load dataset_tecator.mat;
lambda = linspace(850, 1050, 100);
Qlabels = {'Moisture', 'Fat', 'Protein'};

m = size(X, 1);
p = size(X, 2);
q = size(Y, 2);
r = min(p,q);

MSE = nan(r, r);
for k=1:r
  for l=1:min(k, r)
    [P, D, Q, muX, muY, E, Fmax, MSEcv, Fmaxcv, A, B, W] = ...
        moxregress(X, Y, k, l, 'CV', 10, 'MCReps', 10);
    MSE(k,l) = MSEcv;
  end;
end;


figure(10);
bar3(MSE);
ylabel('{\itk}');
xlabel('{\itl}');
zlabel('MSE (cross-validated)');


k = r;
l = 2;
[P, D, Q, muX, muY, E, Fmax, A, B, W] = mox(X, Y, k, l);

fig1 = figure(1);
%[axs, lgd] = moxplot(P, D, Q, X - ones(m, 1)*muX, Y - ones(m, 1)*muY, 'TECATOR', [], [], Qlabels, [], [], [], true, lambda);
[axs, lgd] = moxplot(P, D, Q, X - ones(m, 1)*muX, Y - ones(m, 1)*muY, ...
    'TECATOR', [], Qlabels, [], [], true, lambda);
xlabel(axs(1), '\lambda [nm]');
xlabel(axs(3), '');
currentPosition = fig1.Position; 
fig1.Position = [currentPosition(1), currentPosition(2), 0.9*currentPosition(3), 1.6*currentPosition(4)];

% Set up printing options
set(fig1, 'Color', 'white');   % Ensure the figure background is white (optional)
set(fig1, 'Renderer', 'painters');
%print(fig1, '-depsc', 'img/mox_tecator.eps');