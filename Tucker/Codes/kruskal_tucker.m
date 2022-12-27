function X = kruskal_tucker(G, A, B, C)
[I, P] = size(A);
[J, Q] = size(B);
[K, R] = size(C);
X = zeros(I, J, K);
for k = 1:K
    k
    for j = 1:J
        for i = 1:I
            for p = 1:P
                for q = 1:Q
                    for r = 1:R
                        X(i, j, k) = X(i, j, k) + G(p, q, r) * A(i, p) * B(j, q) * C(k, r);
                    end
                end
            end
        end
    end
end