%%%Code to generate the comparison of RPCA algorithms for simulated data


clear;
clc;
close all

addpath('YALL1_v1.4')

% pctrunonall warning_off_etc


%% Data Generationpctrunonall warning_off_etc

n = 1000;
t_max = 12000;
t_train = 100;
miss_s = 0;
miss_s_pca = 0;
alpha = 300;
f = 100;
MC = 1;

t_calc_pca = [alpha - 1 : alpha : t_max - t_train];

sigma = 1e-5;
theta_degree = 15;

t_reprocs_off = 0;
%t_ncrpca = 0;

outfracrowrange = [0.01, 0.3, 0.8];
rrange = [3, 10, 30];

% PhaseTransNORST = zeros(length(outfracrowrange), length(rrange));
% PhaseTransAltProj = zeros(length(outfracrowrange), length(rrange));

PhaseTransNORST = zeros(length(outfracrowrange), length(rrange), MC);
PhaseTransAltProj = zeros(length(outfracrowrange), length(rrange), MC);



ttall = tic;

for bb  = 1 : length(outfracrowrange)
    b0 = outfracrowrange(bb);
    
    for rr = 1 : length(rrange)
        r_0 = rrange(rr);
        
        %         temp_err_L_off = zeros(MC, t_max - t_train);
        %         temp_err_L_ncrpca = zeros(MC, t_max- t_train);
        %
        %         temp_SE_reprocs_off = zeros(MC, ceil((t_max- t_train)/alpha) - 1);
        %         temp_SE_ncrpca = zeros(MC, ceil((t_max- t_train)/alpha) - 1);
        parfor mc = 1 : MC
            fprintf('Monte-Carlo iteration %d\t b0 %.1f\t r %d', mc, b0, r_0);
            
            %%%Generating support set and sparse vectors
            T = zeros(n, t_max);
            rho_train = 0.02;
            rho = b0;
            x_min = 10;
            x_max = 20;
            
            BernMat = rand(n, t_max);
            T(:, 1 : t_train) = 1.* (BernMat(:, 1 : t_train) <= rho_train);
            T(:, t_train + 1 : end) = 1 .* (BernMat(:, t_train + 1 : t_max) <= rho);
            S = (x_min + (x_max - x_min) * rand(n, t_max)) .* T;
            
            
            
            %%%Generate low-rank matrix
            %r_0 = 30;
            r_1 = 0;
            r_2 = 0;
            r = r_0;
            L = zeros(n, t_max);
            %     diag_entries = [linspace(sqrt(f), sqrt(f)/2, r_0 - r_1), ...
            %         ones(1 , r_1)];
            
            diag_entries = [linspace(sqrt(f), sqrt(f)/2, r_0)];
            t_1 = 3000;
            t_2 = 8000;
            P = orth(randn(n, r_0));
            coeff_train = zeros(r_0, t_max);
            
            for cc = 1 : r_0
                coeff_train(cc, :) = -diag_entries(cc) + ...
                    2 * diag_entries(cc) * rand(1, t_max);
            end
            
            Btemp1 = randn(n);
            B1 = (Btemp1 - Btemp1')/2;
            Btemp2 = randn(n);
            B2 = (Btemp2 - Btemp2')/2;
            
            delta1 = .5e-3;
            delta2 = 0.8 * delta1;
            
            PP1 = expm(delta1 * B1)  * P;
            PP2 = expm(delta2 * B2) * PP1;
            
            L(:, 1:t_1) = P(:, 1:r_0) * coeff_train(:, 1:t_1);
            L(:, t_1+1:t_2) = PP1 * coeff_train(:, t_1+1:t_2);
            L(:, t_2 + 1 : end) = PP2 * coeff_train(:, t_2+1:end);
            M = L + S;
            
            %% Main online robust PCA algorithm section
            
            %%%Algorithm parameters
            K = 8;
            omega = x_min / 2;
            ev_thresh = 7.5961e-04;
            
            %%%Call to offline NORST
            fprintf('\tOffline NORST\t');
            %             tt7 = tic;
            P_init = orth(ncrpca(M(:, 1 : t_train), r_0, 1e-2, 15));
            [L_hat_off, P_hat_off, S_hat_off, T_hat_off, t_hat_off, ...
                P_track_full_off, P_track_new_off] = ...
                Offline_NORST(M(:, t_train + 1 :end),...
                P_init, ev_thresh, alpha, K, omega);
            %             t_reprocs_off = t_reprocs_off + toc(tt7);
            
            %%Call to AltProj
            fprintf('AltProj\n');
            [L_hat_ncrpca, S_hat_ncrpca] = ncrpca(M(:, t_train + 1 : end), 3 * r_0, 1e-6, 100);
            
            %%Compute performance metrics
            
            temp_SE_reprocs_off = zeros(1, length(t_calc_pca));
            temp_SE_ncrpca = zeros(1, length(t_calc_pca));
            for jj = 1 : length(t_calc_pca)
                P_hat_ncrpca = orth(L_hat_ncrpca(:, 1 : jj));
                if (t_calc_pca(jj) +t_train < t_1)
                    temp_SE_reprocs_off(jj) = ...
                        Calc_SubspaceError(P_track_full_off{1}, ...
                        P(:, 1:r_0));
                    temp_SE_ncrpca(jj) = ...
                        Calc_SubspaceError(P_hat_ncrpca, ...
                        P(:, 1 :r_0));
                elseif((t_calc_pca(jj) +t_train >= t_1) && (t_calc_pca(jj) + t_train < t_2))
                    temp_SE_reprocs_off(jj) = ...
                        Calc_SubspaceError(P_track_full_off{2}, PP1);
                    temp_SE_ncrpca(jj) = ...
                        Calc_SubspaceError(P_hat_ncrpca, PP1);
                else
                    temp_SE_reprocs_off(jj) = ...
                        Calc_SubspaceError(P_track_full_off{3}, PP2);
                    temp_SE_ncrpca(jj) = ...
                        Calc_SubspaceError(P_hat_ncrpca, PP2);
                end
            end
            
%             temp_SE_NORST(mc) = mean(temp_SE_reprocs_off);
%             temp1_SE_ncrpca(mc) = mean(temp_SE_ncrpca);
            
            PhaseTransNORST(bb, rr, mc) = mean(temp_SE_reprocs_off);
            PhaseTransAltProj(bb, rr, mc) = mean(temp_SE_ncrpca);
            
            
        end
        
        %         err_L_off = mean(temp_err_L_off, 1);
        %err_L_ncrpca = mean(temp_err_L_ncrpca, 1);
        
        %         err_SE_reprocs_off = mean(temp_SE_reprocs_off, 1);
        %err_SE_ncrpca = mean(temp_SE_ncrpca, 1);
        
%         PhaseTransNORST(bb, rr) = mean(temp_SE_NORST);
%         PhaseTransAltProj(bb, rr) = mean(temp1_SE_ncrpca);
    end
end
fprintf('\n')
toc(ttall)



%t_reprocs_off = t_reprocs_off / (MC * t_max)
% t_ncrpca = t_ncrpca / (MC * t_max)


% figure;
% subplot(311);
% imagesc(S);
% subplot(312);
% imagesc(S_hat_pca);
% subplot(313);
% imagesc(S_hat_off);

% stry ='$$SE(\hat{P}, P)$$';
% strx = '$$t$$';
% figure
% plot(t_calc_pca, log10(err_SE_reprocs_off), 'gs--', 'LineWidth', 2)
% % plot(t_calc_pca, log10(temp_SE_reprocs_pca(3, :)), 'ks--', 'LineWidth', 2)
% % plot(t_calc_pca, log10(temp_SE_reprocs_off(1, :)), 'b^-', 'LineWidth', 2)
% % plot(t_calc_pca, log10(temp_SE_reprocs_off(2, :)), 'g^-', 'LineWidth', 2)
% % plot(t_calc_pca, log10(temp_SE_reprocs_off(3, :)), 'k^-', 'LineWidth', 2)
%
% % plot(t_calc, log10(err_SE_ncrpca), 'g', 'LineWidth', 2)
% xlabel(strx, 'Interpreter', 'LaTeX', 'FontSize', 14);
% ylabel(stry, 'Interpreter', 'LaTeX', 'FontSize', 14);
% axis tight
%legend('NORST x_{min} = 0.5', 'NORST x_{min} = 5', 'NORST x_{min} = 10', ...
%'Offline x_{min} = 0.5', 'Offline x_{min} = 5', 'Offline x_{min} = 10');
%title('NORST')


% stry ='$$SE(\hat{P}, P)$$';
% strx = '$$t$$';
% figure
% plot(t_calc_pca, log10(temp_SE_reprocs_off(1, :)), 'bs--', 'LineWidth', 2)
% hold
% plot(t_calc_pca, log10(temp_SE_reprocs_off(2, :)), 'gs--', 'LineWidth', 2)
% plot(t_calc_pca, log10(temp_SE_reprocs_off(3, :)), 'ks--', 'LineWidth', 2)
% % plot(t_calc, log10(err_SE_ncrpca), 'g', 'LineWidth', 2)
% xlabel(strx, 'Interpreter', 'LaTeX', 'FontSize', 14);
% ylabel(stry, 'Interpreter', 'LaTeX', 'FontSize', 14);
% axis tight
% legend('x_min = 1', 'x_min = 5', 'x_min = 10');
% title('Offline NORST')
