%%%Demo to implement the NORST-ReProCS algorithm for simulated data.
%This folder contains the code accompanying paper. Please cite this if you
%use the code
%
%[1] "Nearly Optimal Robust Subspace Tracking", Praneeth Narayanamurthy and Namrata Vaswani, ICML, 2018.
%

clear;
clc;
close all
addpath('YALL1_v1.4/')

tic
%% Data Generation
n = 1000;
t_max = 5000;
s = ceil(0.05 * n);
t_train = 200;
miss_s = 0;
alpha = 200;
alpha1 = 100;
f = 1;
cnt = 1;
MC = 1;
err_t = zeros(MC, 16);

%% varing initial error and angles
sigmarange = [1e-3, 1e-5, 1e-7, 1e-9];
for ss = 1 : length(sigmarange)
    sigma = sigmarange(ss);
    fprintf('log(Sigma): %d \n', log10(sigma));
    temp_err_L = zeros(MC, t_max - t_train);
    temp_err_SE = zeros(MC, ceil((t_max - t_train)/alpha));
    
    for mc = 1 : MC
        
        fprintf('Monte-Carlo iteration %d in progress\n', mc);
        
        %%%Generating support set and sparse vectors
        S = zeros(n, t_max);
        rho = 1;
        b0 = 0.1;
        beta = ceil(b0 * alpha1);
        x_max = 25;
        x_min = 10;
        alpha1 = 100;
        num_changes = floor((t_max -t_train)/beta);
        
        num_changes1 = floor(alpha1 / beta);
        
        flag = 0;
        ii1 = 1;
        fval1 = 0;
        for ii = 1 : num_changes
            if(~flag)   %%downward motion
                if(ii1 <= num_changes1)
                    bind = fval1 + (ii1 - 1) * s/rho + 1;
                    sind = min(bind - 1 + s, n);
                    ii1 = ii1 + 1;
                    if(ii1 == num_changes1 + 1)
                        flag = 1;
                        ii1 = 1;
                        fval2 = bind;
                    end
                end
            else
                if(ii1 <= num_changes1)
                    bind = max(fval2 - (ii1 - 1) * s/rho , 1);
                    sind = bind - 1 + s;
                    ii1 = ii1 + 1;
                    if(ii1 == num_changes1 + 1)
                        flag = 0;
                        ii1 = 1;
                    end
                end
            end
            idx = bind : sind;
            jdx = t_train + (ii-1) * beta + 1 : t_train + ii * beta;
            S(idx, jdx) = x_min + ...
                (x_max - x_min) * rand(length(idx), beta);
            T(idx, jdx) = 1;
        end
        
        %%%Generate low-rank matrix
        r = 30;
        
        P = orth(randn(n, r));
        gamma = 1e-3;
        
        Btemp = randn(n);
        B = Btemp - Btemp';
        
        PP1 = expm(gamma * B) * P;
        
        
        L = zeros(n, t_max);
        %             diag_entries = [linspace(sqrt(f), sqrt(f)/2, r_0 - r_1), ...
        %                 ones(1 , r_1)];
        diag_entries = [sqrt(f) * ones(1, r), 1];
        
        
        t_1 = 2500;
        coeff_train = zeros(r, t_max);
        
        for cc = 1 : r
            coeff_train(cc, :) = -diag_entries(cc) + ...
                2 * diag_entries(cc) * rand(1, t_max);
        end
        
        L(:, 1:t_1) = P * coeff_train(:, 1:t_1);
        L(:, t_1 + 1 : end) = PP1 * coeff_train(:, t_1+1:end);
        M = L + S;
        
        %% Main online robust PCA algorithm section
        
        %%%Algorithm parameters
        K = 5;
        omega = x_min/2;
        %jugaad for init errors
        P_init = orth(P + sigma * randn(n, r)/2.5);
        ev_thresh = 7.5961e-04;
        %             ev_thresh = 0.1/3 * (sin(theta1))^2;
        
        %%%Call to online RPCA function
        [L_hat, P_hat, S_hat, T_hat, t_hat, P_track_full, t_calc] = ...
            NORST(M(:, t_train + 1 : end), P_init, ...
            ev_thresh, alpha, K, omega);
        
        %%Compute performance metrics
        temp_err_L(mc, :) = ...
            sqrt(mean((L(:, t_train + 1 : end) - L_hat).^2, 1)) ./ ...
            sqrt(mean(L(:, t_train + 1 : end).^2, 1));
        miss_s = ...
            miss_s + (length(find(S_hat))- length(find(S)))/numel(S);
        
        %%Calculate the subspace error
        for jj = 1 : length(t_calc)
            if (t_calc(jj) + t_train < t_1)
                temp_SE_Phat_P(mc, jj) = ...
                    Calc_SubspaceError(P_track_full{jj}, P);
            elseif (t_calc(jj) + t_train >= t_1)
                temp_SE_Phat_P(mc, jj) = ...
                    Calc_SubspaceError(P_track_full{jj}, PP1);
            end
        end
        %fprintf('\n\n');
    end
    err_L(cnt, :) = mean(temp_err_L, 1);
    SE_Phat_P(cnt, :) = mean(temp_SE_Phat_P, 1);
    cnt = cnt + 1;
end
toc

%%call to this works only when there are 16 sigma-theta combination, otherwise need to manually generate figures
%FigGenReProCS

