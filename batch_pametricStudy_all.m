clc; clear; close all;

% Latent variable information
d_XY   = 2;                  % co-varying dimensions shared by X and Y
eps_y  = [8, 4, 2, 1, 0];
eps_x  = [8, 4, 2, 1, 0];
d_X    = eps_x*d_XY;         % co-varying dimensions in X only
d_Y    = eps_y*d_XY;         % co-varying dimensions in Y only

% Number of observations, predictors, and responses
m_vec = [1, 2, 4]*40;   % no. of observations
q_vec = [20, 34];       % d_X+d_XY; % features in responses (multivariate output array)
p_vec = [1, 2, 4]*20;   % features in predictor (multivariate input array)
w_amp = 0.5;            % noise level

% Reproducibility parameters: Repitions, MCreps, CV
n_random_seeds = 20;  % Random numbers
n_repetitions  = 50;  % Number of Monte Carlo repetitions over CV
CV             = 10;  % 10-folsd cross-validation

for im=1:length(m_vec)
    m = m_vec(im);
    for iq=1:length(q_vec)
        q=q_vec(iq);

        for jp=1:length(p_vec)
            p=p_vec(jp);
            filename   = sprintf('parametric_study_m%d_q%d_p%d.xlsx', m, q, p);
            %fname_nsr = sprintf('NSR_values_m%d_q%d_p%d.xlsx', m, q, p);

            for iey = 1:length(d_Y)
                for iex = 1:length(d_X)
                    d      = d_X(iex) + d_XY + d_Y(iey); % total dimensionality of fluctuations
                    idx_X  = (1:d_X(iex));
                    idx_XY = d_X(iex) + (1:d_XY);
                    idx_Y  = d_X(iex) + d_XY + (1:d_Y(iey));

                    % Compute MSE for all
                    benchmark_vishal;

                    %--- Write CV MSE over (m,p,q,ex,ey) at dXY, dX+dXY, and at minGlobal
                    MSE_all = MSEvalues_all(d_X(iex), d_XY, d_Y(iey), ...
                        MSEmox_mean, MSEmox_kisl_mean, ...
                        MSEpls_mean, MSEcca_mean, MSEols_mean);
                    label   = sprintf('(%d, %d)', eps_x(iex), eps_y(iey));
                    rowData = [{label}, num2cell(MSE_all)];
                    writecell(rowData, filename, 'WriteMode', 'append');

                    % Write NSR values over (m,p,q,ex,ey) at dXY
                    %[P,D,Q,~,~, ~, ~, ~, ~] = mox(X, Y, k, d_XY);
                    %[~, ~, ~, ~, BETA]      = plsregress(X, Y, d_XY);
                    % % MOX_NSR
                    %Ypred_mox      = X*P*D*Q';
                    %Residual_mox   = Y - Ypred_mox;
                    %noise_var_mox  = var(Residual_mox(:));
                    %signal_var_mox = var(Ypred_mox(:));
                    %NSR_mox        = noise_var_mox/signal_var_mox;
                    % % PLS_NSR
                    %Ypred_pls      = [ones(size(X, 1), 1), X] * BETA;
                    %Residual_pls   = Y - Ypred_pls;
                    %noise_var_pls  = var(Residual_pls(:));
                    %signal_var_pls = var(Ypred_pls(:));
                    %NSR_pls        = noise_var_pls/signal_var_pls;
                    % % Write
                    %NSR_all        = [NSR_dY, NSR_mox, NSR_pls];
                    %rowData_nsr    = [{label}, num2cell(NSR_all)];
                    %writecell(rowData_nsr, fname_nsr, 'WriteMode', 'append');

                    %--- Print messages
                    fprintf('\n (iex,iey) = (%d, %d)', eps_x(iex), eps_y(iey));
                    %fprintf('.       ');
                    %fprintf('Noise-to-Signal Ratio ...(NSR_dY, NSR_mox, NSR_pls) = (%.2f, %.2f, %.2f) \n', NSR_dY, NSR_mox, NSR_pls);
                    fprintf('\n');
                end % eps_x
                % Write blank row
                blankRow = num2cell(NaN(1, size(rowData, 2)));
                writecell(blankRow, filename, 'WriteMode', 'append');
                %blankRow_nsr = num2cell(NaN(1, size(rowData_nsr, 2)));
                %writecell(blankRow_nsr, fname_nsr, 'WriteMode', 'append');
            end  % eps_y
        end      % p_vec
        fprintf('(p,q)=(%d,%d) \n', p, q)
    end          % q_vec
    fprintf('Done: m = %d \n', m)
    disp('----------------------------------------------')
    disp('')
    disp('')
end % m_vec