clc; clear all; close all;

% Latent variable information
d_XY   = 2;                  % co-varying dimensions shared by X and Y
eps_x  = [8, 4, 2, 1, 0];
eps_y  = [8, 4, 2, 1, 0];
d_X    = eps_x*d_XY;         % co-varying dimensions in X only
d_Y    = eps_y*d_XY;         % co-varying dimensions in Y only

% based on k<=min(p,q). q is set according to it. For p two cases [p>p, p=q].
% The third case: p<q will violate k<min(p,q).
q_vec = [20];         % d_X+d_XY; % features in responses (multivariate output array)
p_vec = [20, 40, 80]; % features in predictor (multivariate input array)
m_vec = [1,2]*40;     % no. of observations
w_amp = 0.5;          % noise level

% standard: [20, 50, 10]
n_random_seeds = 20;
n_repetitions  = 50;  % Number of Monte Carlo repetitions over CV
CV             = 10;  % 10-folsd cross-validation

for im=1:length(m_vec)
    m = m_vec(im);

    for iq=1:length(q_vec)
        q=q_vec(iq);

        for jp=1:length(p_vec)
            %p=p_vec(jp)*q;    % use when p=>q
            p=p_vec(jp);%*q;   % use when p<q
            
						filename = sprintf('parametric_study_m%d_q%d_p%d.xlsx', m, q, p);

            %for iex = iq:length(d_X)
            for iey = 1:length(d_Y)
                %d_Y = d_Y_vec(iey);

                %for iey = 1:length(d_Y)
								for iex = 1:length(d_X)
                    %d_X = d_X_vec(iex);

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
                    fprintf('(iex,iey) = (%d, %d) \n', eps_x(iex), eps_y(iey))
                    fprintf('\n')
                end % eps_x
            end     % eps_y
        end         % p_vec
				fprintf('(p,q)=(%d,%d) \n', p, q)
    end             % q_vec
		fprintf('m = %d \n', m)
		disp('----------------------------------------------')
		disp('')
		disp('')
end                 % m_vec
