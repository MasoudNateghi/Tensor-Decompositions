function t_norm = tensor_norm(T)
I = size(T, 1);
J = size(T, 2);
K = size(T, 3);
t_norm = 0;
for k = 1:K
    for j = 1:J
        for i = 1:I
            t_norm = t_norm + T(i, j, k) ^ 2;
        end
    end
end
t_norm = sqrt(t_norm);