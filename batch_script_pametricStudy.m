clc;

% Simulation parameters
m              = 80;  % no. of observations
p              = 60;  % features in predictor (multivariate input array)
q              = 20;  % features in responses (multivariate output array)
w_amp          = 0.5; % noise level
% standard: [20, 50, 10]
n_random_seeds = 20;
n_repetitions  = 50;  % Number of Monte Carlo repetitions over CV
CV             = 10;  % 10-folsd cross-validation

% Filename
% filename = 'parametric_study.xlsx';
filename = sprintf('parametric_study_m%d_p%d_q%d.xlsx', m, p, q);

% Latent variable information
d_XY   = 2;                  % co-varying dimensions shared by X and Y
eps_x  = flip([0, 1, 2, 4, 8]);
eps_y  = flip([0, 1, 2, 4, 8]);
d_X    = eps_x*d_XY;         % co-varying dimensions in X only
d_Y    = eps_y*d_XY;         % co-varying dimensions in Y only

for iey = 1:length(eps_y)
    for iex = 1:length(eps_x)

        % compute idx_[]
        d      = d_X(iex) + d_XY + d_Y(iey); % total dimensionality of fluctuations
        idx_X  = (1:d_X(iex));
        idx_XY = d_X(iex) + (1:d_XY);
        idx_Y  = d_X(iex) + d_XY + (1:d_Y(iey));

        % Compute MSE for all
        benchmark_vishal;

        % Write varying parameters and MSE levels
        MSE_all = MSEvalues_all(d_X(iex), d_XY, d_Y(iey), MSEmox_mean, MSEmox_kisl_mean, ...
            MSEpls_mean, MSEcca_mean, MSEols_mean);
        label   = sprintf('(%d, %d)', eps_x(iex), eps_y(iey));
        rowData = [{label}, num2cell(MSE_all)];
        writecell(rowData, filename, 'WriteMode', 'append');
				% black row
        blankRow = repmat({''}, 1, size(rowData, 2));
        writecell(blankRow, filename, 'WriteMode', 'append');
        fprintf('completed (iex,iey) = (%d, %d) \n', iex, iey)
        fprintf('\n \n')
    end
    writecell(repmat({''}, 1, size(rowData,2)), filename, 'WriteMode', 'append');
end
% add two more blank rows at the end of all (eps_x, eps_y) combination block
writecell(repmat({''}, 1, size(rowData,2)), filename, 'WriteMode', 'append');
writecell(repmat({''}, 1, size(rowData,2)), filename, 'WriteMode', 'append');
writecell(repmat({''}, 1, size(rowData,2)), filename, 'WriteMode', 'append');
