clear; close all; clc;
X = zeros(112, 92, 50);
cd ORL
count = 1;
dirs = dir;
for i = 3:length(dirs)
    s = convertCharsToStrings(dirs(i).name);
    for k = 1:10
        X(:, :, count) = im2double(imread(s+"\"+num2str(k)+".pgm"));
        count = count + 1;
    end
end
cd ..
%% als (developed)
% !!!!!!!!!!!!! LONG RUNTIME !!!!!!!!!!!!!
% [G, A1, A2, A3] = hooi(X, 5, 5, 5);
for i = 1:5
    plot(A3(:, i))
    hold on
end
xlabel("People", "Interpreter","latex")
legend(["s1", "s2", "s3", "s4", "s5"], "interpreter", "latex")
title("HOOI (developed) in 10th iteration!!", "Interpreter","latex")
grid on
%% tucker_als
T = tucker_als(tensor(X), [5, 5, 5]);
U3 = T.U{3};
for i = 1:5
    plot(U3(:, i))
    hold on
end
grid on
xlabel("People", "Interpreter","latex")
legend("s1", "s2", "s3", "s4", "s5", "interpreter", "latex")
title("HOOI (Toolbox) in 10th, R = [5, 5, 5]", "Interpreter","latex")
%% ntd
opts.maxit = 500;
[A,C] = ntd(tensor(X), [35, 35, 5], opts);
U3 = A{3};
for i = 1:5
    plot(U3(:, i))
    hold on
end
grid on
xlabel("People", "Interpreter","latex")
legend("s1", "s2", "s3", "s4", "s5", "interpreter", "latex")
title("Non-negative Tensor Decomposition, R = [35, 35, 5]", "Interpreter","latex")




















