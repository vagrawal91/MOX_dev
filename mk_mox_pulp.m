%load ../../modes/analysis/modes_kvarnsveden64.mat
load('/Users/vishalagrawal/DriveD/mox_sim/sim_dimsame/data/milldata_kvarnsveden64.mat')
close all;

%idxX = [1 2 3 4 5 6 7 9 10 12 13 15 16]; % 8 11 14 and 17 excluded
idxX = [1 5 6 7 8 9 10 11 12 13 14 15 16 17];
idxY = [1 2 3 4 5 6 7 8]; % 9 excluded
X = X2z(:,idxX);
Y = X3z(:,idxY);
m = size(X, 1);
p = size(X, 2);
q = size(Y, 2);
r = min(p,q);

MSE = nan(r, r);
for k=1:r
  for l=1:min(k, r)
    [P, D, Q, muX, muY, E, Fmax, MSEcv, Fmaxcv, A, B, W] = moxregress(X, Y, k, l, 'CV', 10, 'MCReps', 10);
    MSE(k,l) = MSEcv / q;
  end
end

figure(10);
bar3(MSE);
ylabel('{\itk}');
xlabel('{\itl}');
zlabel('MSE (cross-validated)');

k = q;
l = 3;
[P, D, Q, muX, muY, E, Fmax, A, B, W] = mox(X, Y, k, l);
flipidx = (Q(3,:) < 0);  % Change sign to get positive TI
P(:,flipidx) = -P(:,flipidx);
Q(:,flipidx) = -Q(:,flipidx);

fig1 = figure(1);
[axs, lgd] = moxplot(P, D, Q, X - ones(m, 1)*muX, Y  - ones(m, 1)*muY, 'Pulp', X2lbl(idxX), X3lbl(idxY));
xlabel(axs(1), '');
xlabel(axs(2), '');

fsize = [12 17]; % cm
set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 fsize], 'PaperSize', fsize);
set(gcf, 'Units', 'centimeters', 'Position', [0 0 fsize]);
%exportgraphics(gcf, 'img/mox_pulp.pdf', 'ContentType', 'vector');