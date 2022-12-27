function [U1, U2, U3] = cp_als(T, U1_0, U2_0, U3_0)
numItter = 50;
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
U1 = U1_0; U2 = U2_0; U3 = U3_0;
for i = 1:numItter
    V1 = (U2' * U2) .* (U3' * U3);
    U1 = T1 * khatri_rao(U3, U2) * pinv(V1);
    V2 = (U1' * U1) .* (U3' * U3);
    U2 = T2 * khatri_rao(U3, U1) * pinv(V2);
    V3 = (U1' * U1) .* (U2' * U2);
    U3 = T3 * khatri_rao(U2, U1) * pinv(V3);
end
end

function AB = khatri_rao(A, B)
R = size(A, 2);
AB = zeros(size(A, 1) * size(B, 1), R);
for i = 1:R
    temp = B(:, i) * A(:, i)';
    AB(:, i) = temp(:);
end
end