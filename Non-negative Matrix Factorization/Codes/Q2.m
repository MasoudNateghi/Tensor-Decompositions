clear; close all; clc;
load("swimmer.mat")
Y = zeros(256, 9*14);
for i = 1:256
    Y(i, :) = reshape(A{1, i}, 1, 9*14);
end
imagesc(Y)
colormap jet
e_ALS = zeros(1, 9*14);
e_MULT = zeros(1, 9*14);
for j = 1:9*14
    [B1, C1] = nnmf(Y, j, "algorithm", "als");
    [B2, C2] = nnmf(Y, j, "algorithm", "mult");
    e_ALS(j) = norm(Y - B1 * C1, "fro");
    e_MULT(j) = norm(Y - B2 * C2, "fro");
end
figure;
plot(1:9*14, e_ALS)
xlabel("j", "Interpreter","latex")
ylabel("$$||E||_F$$", "Interpreter","latex")
title("Errors Using ALS algorithm", "Interpreter","latex")
figure;
plot(1:9*14, e_MULT)
xlabel("j", "Interpreter","latex")
ylabel("$$||E||_F$$", "Interpreter","latex")
title("Errors Using Multiplicative algorithm", "Interpreter","latex")
%%
j_best = 16;
[B, C] = nnmf(Y, j_best, "algorithm","mult");
figure;
for i = 1:j_best
    subplot(4, 4, i)
    imagesc(reshape(C(i, :), 9, 14))
end
colormap jet
figure;
imagesc(B)
colormap jet
title("B", "Interpreter","latex")
figure;
imagesc(Y - B * C)
colormap jet
title("E", "Interpreter","latex")