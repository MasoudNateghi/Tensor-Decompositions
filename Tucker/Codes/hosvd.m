function [U1_0, U2_0, U3_0] = hosvd(T, R1, R2, R3)
% unfolding mode 1
T1 = unfold(T, 1);

% unfolding mode 2
T2 = unfold(T, 2);

% unfolding mode 3
T3 = unfold(T, 3);

[U, ~, ~] = svd(T1);
U1_0 = U(:, 1:R1);

[U, ~, ~] = svd(T2);
U2_0 = U(:, 1:R2);

[U, ~, ~] = svd(T3);
U3_0 = U(:, 1:R3);