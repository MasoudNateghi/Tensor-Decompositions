clear; close all; clc;
I = 6; 
J = 4;
K = 2;
R = 5;
U1_org = randn(I, R);
U2_org = randn(J, R);
U3_org = randn(K, R);
T = zeros(I, J, K);
% construct tensor from outer product
for k = 1:K
    for j = 1:J
        for i = 1:I
            for r = 1:R
                T(i, j, k) = T(i, j, k) + U1_org(i, r) * U2_org(j, r) * U3_org(k, r);
            end
        end
    end
end
U1_0 = randn(I, R);
U2_0 = randn(J, R);
U3_0 = randn(K, R);
[U1, U2, U3] = cp_als(T, U1_0, U2_0, U3_0);
T_hat = zeros(I, J, K);
for k = 1:K
    for j = 1:J
        for i = 1:I
            for r = 1:R
                T_hat(i, j, k) = T_hat(i, j, k) + U1(i, r) * U2(j, r) * U3(k, r);
            end
        end
    end
end