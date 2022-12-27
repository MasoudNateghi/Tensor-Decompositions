%% part a/b
clear; close all; clc;
SNR = 0:10:60;
error_cp_als = zeros(length(SNR), 1);
error_cpd_als = zeros(length(SNR), 1);
error_cpd_sd = zeros(length(SNR), 1);
error_cpd_minf = zeros(length(SNR), 1);

error_cp_als_HOSVD = zeros(length(SNR), 1);
error_cpd_als_HOSVD = zeros(length(SNR), 1);
error_cpd_sd_HOSVD = zeros(length(SNR), 1);
error_cpd_minf_HOSVD = zeros(length(SNR), 1);

for m = 1:50
    m %#ok<NOPTS> 

    % principal loading factors
    U1_org = randn(6, 3);
    U2_org = randn(6, 3);
    U3_org = randn(6, 3);
    U_org = {U1_org, U2_org, U3_org};
    T = zeros(6, 6, 6);

    % construct tensor from outer product
    for k = 1:6
        for j = 1:6
            for i = 1:6
                for r = 1:3
                    T(i, j, k) = T(i, j, k) + U1_org(i, r) * U2_org(j, r) * U3_org(k, r);
                end
            end
        end
    end

    % noise addition
    N = randn(6, 6, 6);
    X = zeros(6, 6, 6, 4);
    count = 1;
    for i = 1:length(SNR)
        alpha = (tensor_norm(T) ^ 2) / (10 ^ (SNR(i)/10)) / (tensor_norm(N) ^ 2);
        X(:, :, :, count) = T + sqrt(alpha) * N;
        count = count + 1;
    end

    % initial loading matrix (random)
    U1_0 = randn(6, 3);
    U2_0 = randn(6, 3);
    U3_0 = randn(6, 3);
    U0 = {U1_0, U2_0, U3_0};

    for i = 1:length(SNR)
        % CP decomposition using developed developed algorithm
        [U1, U2, U3] = cp_als(X(:, :, :, i), U1_0, U2_0, U3_0);
        esU = {U1, U2, U3};
        error_cp_als(i) = error_cp_als(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd_als(X(:, :, :, i), U0);
        error_cpd_als(i) = error_cpd_als(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd3_sd(X(:, :, :, i), U0);
        error_cpd_sd(i) = error_cpd_sd(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd_minf(X(:, :, :, i), U0);
        error_cpd_minf(i) = error_cpd_minf(i) + TMSFE(U_org, esU) / 50;

        % initial loading matrix (HOSVD)
        [U1_0_HOSVD, U2_0_HOSVD, U3_0_HOSVD] = hosvd(X(:, :, :, i), 3);
        U0_HOSVD = {U1_0_HOSVD, U2_0_HOSVD, U3_0_HOSVD};
        
        % CP decomposition using developed developed algorithm
        [U1, U2, U3] = cp_als(X(:, :, :, i), U1_0_HOSVD, U2_0_HOSVD, U3_0_HOSVD);
        esU = {U1, U2, U3};
        error_cp_als_HOSVD(i) = error_cp_als_HOSVD(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd_als(X(:, :, :, i), U0_HOSVD);
        error_cpd_als_HOSVD(i) = error_cpd_als_HOSVD(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd3_sd(X(:, :, :, i), U0_HOSVD);
        error_cpd_sd_HOSVD(i) = error_cpd_sd_HOSVD(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd_minf(X(:, :, :, i), U0_HOSVD);
        error_cpd_minf_HOSVD(i) = error_cpd_minf_HOSVD(i) + TMSFE(U_org, esU) / 50;
    end
end

figure;
plot(0:10:60, error_cp_als)
hold on
plot(0:10:60, error_cpd_als)
plot(0:10:60, error_cpd_sd)
plot(0:10:60, error_cpd_minf)
legend("cp ALS (developed)", "cp ALS", "cp SD", "cp MINF", "interpreter", "latex")
xlabel("SNR")
ylabel("TMSFE (Error)")
title("Canonical Polyadic Tensor Decomposition with Random Initialization", "Interpreter","latex")
grid minor

figure;
plot(0:10:60, error_cp_als_HOSVD)
hold on
plot(0:10:60, error_cpd_als_HOSVD)
plot(0:10:60, error_cpd_sd_HOSVD)
plot(0:10:60, error_cpd_minf_HOSVD)
legend("cp ALS (developed)", "cp ALS", "cp SD", "cp MINF", "interpreter", "latex")
xlabel("SNR")
ylabel("TMSFE (Error)")
title("Canonical Polyadic Tensor Decomposition with HOSVD Initialization", "Interpreter","latex")
grid minor

%% part c/d
clear; close all; clc;
SNR = 0:10:60;
error_cp_als = zeros(length(SNR), 1);
error_cpd_als = zeros(length(SNR), 1);
error_cpd_sd = zeros(length(SNR), 1);
error_cpd_minf = zeros(length(SNR), 1);

error_cp_als_HOSVD = zeros(length(SNR), 1);
error_cpd_als_HOSVD = zeros(length(SNR), 1);
error_cpd_sd_HOSVD = zeros(length(SNR), 1);
error_cpd_minf_HOSVD = zeros(length(SNR), 1);

for m = 1:50
    m
    % principal loading factors
    U1_org = zeros(6, 3);
    U1_org(:, 1) = randn(6, 1);
    U1_org(:, 2) = U1_org(:, 1) + 0.5 * randn(6, 1);
    U1_org(:, 3) = randn(6, 1);

    U2_org = zeros(6, 3);
    U2_org(:, 1) = randn(6, 1);
    U2_org(:, 2) = U2_org(:, 1) + 0.5 * randn(6, 1);
    U2_org(:, 3) = randn(6, 1);

    U3_org = randn(6, 3);

    U_org = {U1_org, U2_org, U3_org};

    T = zeros(6, 6, 6);
    
    % construct tensor from outer product
    for k = 1:6
        for j = 1:6
            for i = 1:6
                for r = 1:3
                    T(i, j, k) = T(i, j, k) + U1_org(i, r) * U2_org(j, r) * U3_org(k, r);
                end
            end
        end
    end
    
    % noise addition
    N = randn(6, 6, 6);
    X = zeros(6, 6, 6, 4);
    count = 1;
    for i = 1:length(SNR)
        alpha = (tensor_norm(T) ^ 2) / (10 ^ (SNR(i)/10)) / (tensor_norm(N) ^ 2);
        X(:, :, :, count) = T + sqrt(alpha) * N;
        count = count + 1;
    end

    % initial loading matrix (random)
    U1_0 = randn(6, 3);
    U2_0 = randn(6, 3);
    U3_0 = randn(6, 3);
    U0 = {U1_0, U2_0, U3_0};

    for i = 1:length(SNR)
        % CP decomposition using developed developed algorithm
        [U1, U2, U3] = cp_als(X(:, :, :, i), U1_0, U2_0, U3_0);
        esU = {U1, U2, U3};
        error_cp_als(i) = error_cp_als(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd_als(X(:, :, :, i), U0);
        error_cpd_als(i) = error_cpd_als(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd3_sd(X(:, :, :, i), U0);
        error_cpd_sd(i) = error_cpd_sd(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd_minf(X(:, :, :, i), U0);
        error_cpd_minf(i) = error_cpd_minf(i) + TMSFE(U_org, esU) / 50;

        % initial loading matrix (HOSVD)
        [U1_0_HOSVD, U2_0_HOSVD, U3_0_HOSVD] = hosvd(X(:, :, :, i), 3);
        U0_HOSVD = {U1_0_HOSVD, U2_0_HOSVD, U3_0_HOSVD};
        
        % CP decomposition using developed developed algorithm
        [U1, U2, U3] = cp_als(X(:, :, :, i), U1_0_HOSVD, U2_0_HOSVD, U3_0_HOSVD);
        esU = {U1, U2, U3};
        error_cp_als_HOSVD(i) = error_cp_als_HOSVD(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd_als(X(:, :, :, i), U0_HOSVD);
        error_cpd_als_HOSVD(i) = error_cpd_als_HOSVD(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd3_sd(X(:, :, :, i), U0_HOSVD);
        error_cpd_sd_HOSVD(i) = error_cpd_sd_HOSVD(i) + TMSFE(U_org, esU) / 50;

        % CP decomposition using developed ALS algorithm (tensorlab)
        esU = cpd_minf(X(:, :, :, i), U0_HOSVD);
        error_cpd_minf_HOSVD(i) = error_cpd_minf_HOSVD(i) + TMSFE(U_org, esU) / 50;
    end
end


figure;
plot(0:10:60, error_cp_als)
hold on
plot(0:10:60, error_cpd_als)
plot(0:10:60, error_cpd_sd)
plot(0:10:60, error_cpd_minf)
legend("cp ALS (developed)", "cp ALS", "cp SD", "cp MINF", "interpreter", "latex")
xlabel("SNR")
ylabel("TMSFE (Error)")
title("Canonical Polyadic Tensor Decomposition with Random Initialization", "Interpreter","latex")
grid minor

figure;
plot(0:10:60, error_cp_als_HOSVD)
hold on
plot(0:10:60, error_cpd_als_HOSVD)
plot(0:10:60, error_cpd_sd_HOSVD)
plot(0:10:60, error_cpd_minf_HOSVD)
legend("cp ALS (developed)", "cp ALS", "cp SD", "cp MINF", "interpreter", "latex")
xlabel("SNR")
ylabel("TMSFE (Error)")
title("Canonical Polyadic Tensor Decomposition with HOSVD Initialization", "Interpreter","latex")
grid minor