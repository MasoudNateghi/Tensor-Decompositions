function T_unfold = unfold(T, mode)
[I, J, K] = size(T);
if mode == 1
    % unfolding mode 1
    T_unfold = zeros(I, J*K);
    count = 1;
    for k = 1:K
        for j = 1:J
            T_unfold(:, count) = T(:, j, k);
            count = count + 1;
        end
    end
elseif mode == 2
    % unfolding mode 2
    T_unfold = zeros(J, I*K);
    count = 1;
    for k = 1:K
        for i = 1:I
            T_unfold(:, count) = T(i, :, k);
            count = count + 1;
        end
    end
elseif mode == 3
    % unfolding mode 3
    T_unfold = zeros(K, I*J);
    count = 1;
    for j = 1:J
        for i = 1:I
            T_unfold(:, count) = T(i, j, :);
            count = count + 1;
        end
    end
end