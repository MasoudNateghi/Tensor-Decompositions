function [G, U1, U2, U3] = hooi(T, R1, R2, R3)
[I, J, K] = size(T);
% hosvd initialization
[U1, U2, U3] = hosvd(T, R1, R2, R3);
numItter = 10;
for i = 1:numItter
    i %#ok<NOPRT> 
    Z1 = kruskal_tucker(T, eye(I), U2', U3');
    [U, ~, ~] = svd(unfold(Z1, 1));
    U1 = U(:, 1:R1);

    Z2 = kruskal_tucker(T, U1', eye(J), U3');
    [U, ~, ~] = svd(unfold(Z2, 2));
    U2 = U(:, 1:R2);

    Z3 = kruskal_tucker(T, U1', U2', eye(K));
    [U, ~, ~] = svd(unfold(Z3, 3));
    U3 = U(:, 1:R3);
end
G = kruskal_tucker(T, U1', U2', U3');