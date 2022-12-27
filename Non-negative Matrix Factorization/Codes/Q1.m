clear; close all; clc;
B = rand(6, 3);
C = rand(3, 4);
A = B * C;
for j = 1:4
    % initial values
    B0 = rand(6, j);
    C0 = rand(j, 4);
    [B1, C1] = mynnmf(A, j, "als", B0, C0);
    [B2, C2] = mynnmf(A, j, "mult", B0, C0);
    [B3, C3] = nnmf(A, j, "algorithm", "als", "w0", B0, "h0", C0);
    [B4, C4] = nnmf(A, j, "algorithm", "mult", "w0", B0, "h0", C0);
    disp("j = " + num2str(j))
    disp("B (ALS mynnmf): ")
    disp(B1)
    disp("C (ALS mynnmf): ")
    disp(C1)
    disp("A_hat (ALS mynnmf):")
    disp(B1*C1)
    disp("Error (ALS mynnmf):")
    disp(norm(A-B1*C1, "fro"))
    disp("B (ALS nnmf): ")
    disp(B3)
    disp("C (ALS nnmf): ")
    disp(C3)
    disp("A_hat (ALS nnmf):")
    disp(B3*C3)
    disp("Error (ALS nnmf):")
    disp(norm(A-B3*C3, "fro"))
    disp("B (Multiplicative mynnmf): ")
    disp(B2)
    disp("C (Multiplicative mynnmf): ")
    disp(C2)
    disp("A_hat (Multiplicative mynnmf):")
    disp(B2*C2)
    disp("Error (Multiplicative mynnmf):")
    disp(norm(A-B2*C2, "fro"))
    disp("B (Multiplicative nnmf): ")
    disp(B4)
    disp("C (Multiplicative nnmf): ")
    disp(C4)
    disp("A_hat (Multiplicative nnmf):")
    disp(B4*C4)
    disp("Error (Multiplicative nnmf):")
    disp(norm(A-B4*C4, "fro"))
    disp("======================================================")
end
%%
clear; close all; clc;
% produce A
B = rand(6, 3);
C = rand(3, 4);
E = rand(6, 4);
A = zeros(6, 4, 5);
count = 1;
for SNR = [-10, 0, 10, 30, 50]
    alpha = norm(B * C, "fro") / (10 ^ (SNR / 20) * norm(E, "fro"));
    A(:, :, count) = B * C + alpha * E;
    count = count + 1;
end

e_ALS_mynnmf = zeros(4, 5, 10);
e_ALS_nnmf = zeros(4, 5, 10);
e_MULT_mynnmf = zeros(4, 5, 10);
e_MULT_nnmf = zeros(4, 5, 10);
for k = 1:10
    for j = 1:4
        % initial values
        B0 = rand(6, j);
        C0 = rand(j, 4);
        for i = 1:5
            [B1, C1] = mynnmf(A(:, :, i), j, "als", B0, C0);
            [B2, C2] = mynnmf(A(:, :, i), j, "mult", B0, C0);
            [B3, C3] = nnmf(A(:, :, i), j, "algorithm", "als", "w0", B0, "h0", C0);
            [B4, C4] = nnmf(A(:, :, i), j, "algorithm", "mult", "w0", B0, "h0", C0);
            e_ALS_mynnmf(j, i, k) = norm(A(:, :, i) - B1 * C1, "fro");
            e_MULT_mynnmf(j, i, k) = norm(A(:, :, i) - B2 * C2, "fro");
            e_ALS_nnmf(j, i, k) = norm(A(:, :, i) - B3 * C3, "fro");
            e_MULT_nnmf(j, i, k) = norm(A(:, :, i) - B4 * C4, "fro");
        end
    end
end
e_ALS_nnmf1 = mean(e_ALS_nnmf, 3);
e_MULT_nnmf1 = mean(e_MULT_nnmf, 3);
e_ALS_mynnmf1 = mean(e_ALS_mynnmf, 3);
e_MULT_mynnmf1 = mean(e_MULT_mynnmf, 3);
SNR = [-10, 0, 10, 30, 50];
for i = 1:2
    if i == 1
        e_als = e_ALS_mynnmf1;
        e_mult = e_MULT_mynnmf1;
        str1 = "ALS Algorithm Using mynnmf";
        str2 = "Multiplicative Algorithm Using mynnmf";
    else
        e_als = e_ALS_nnmf1;
        e_mult = e_MULT_nnmf1;
        str1 = "ALS Algorithm Using nnmf";
        str2 = "Multiplicative Algorithm Using nnmf";
    end
    figure;
    for j = 1:4
        plot(SNR, e_als(j, :))
        hold on
    end
    legend("j = 1", "j = 2", "j = 3", "j = 4", "interpreter", "latex")
    xlabel("SNR(dB)", "Interpreter","latex")
    ylabel("$$||E||_F$$", "Interpreter","latex")
    title(str1, "Interpreter","latex")
    figure;
    for j = 1:4
        plot(SNR, e_mult(j, :))
        hold on
    end
    legend("j = 1", "j = 2", "j = 3", "j = 4", "interpreter", "latex")
    xlabel("SNR(dB)", "Interpreter","latex")
    ylabel("$$||E||_F$$", "Interpreter","latex")
    title(str2, "Interpreter","latex")
end