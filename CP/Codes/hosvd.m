function [U1_0, U2_0, U3_0] = hosvd(T, R)
I = size(T, 1);
J = size(T, 2);
K = size(T, 3);
% mode-1 unfolding
T1 = zeros(I, J*K);
count = 1;
for k = 1:K
    for j = 1:J
        T1(:, count) = T(:, j, k);
        count = count + 1;
    end
end
% mode-2 unfolding
T2 = zeros(J, I*K);
count = 1;
for k = 1:K
    for i = 1:I
        T2(:, count) = T(i, :, k);
        count = count + 1;
    end
end
% mode-3 unfolding
T3 = zeros(K, I*J);
count = 1;
for j = 1:J
    for i = 1:I
        T3(:, count) = T(i, j, :);
        count = count + 1;
    end
end
[U1, ~, ~] = svd(T1);
[U2, ~, ~] = svd(T2);
[U3, ~, ~] = svd(T3);
U1_0 = U1(:, 1:R);
U2_0 = U2(:, 1:R);
U3_0 = U3(:, 1:R);