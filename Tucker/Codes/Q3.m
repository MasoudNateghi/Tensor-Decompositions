%% 
% !!!!!!!!!!!!!!!!!! NOTE !!!!!!!!!!!!!!!!!
% Following parts of the code which are related to Tucker Decomposition are
% done with original images (without downsampling). So the runtime is overwhelming. 
%% read data
clear; close all; clc;
X = zeros(10, 9, 32256);
for i = 1:10
    count = 1;
    if i ~= 10
        subject = "yaleB0" + num2str(i);
    else
        subject = "yaleB" + num2str(i);
    end
    for k = ["10", "25", "50", "70"]
        temp = imread("Illumination_Yale\" + subject + "\" + subject + "_P00A-0" + k + "E+00.pgm");
        X(i, count, :) = im2double(temp(:));
        count = count + 1;
    end
    for k = ["10", "25", "50", "70", "00"]
        temp = imread("Illumination_Yale\" + subject + "\" + subject + "_P00A+0" + k + "E+00.pgm");
        X(i, count, :) = im2double(temp(:));
        count = count + 1;
    end
end
%% Tucker Decomposition (part a)
% Patience is the key of success! :))) 

% decomposition
X = tensor(X);
opts.maxit = 5000;
[A1,C1] = ntd(X, [10, 1, 90], opts);
[A3,C3] = ntd(X, [10, 3, 90], opts);
[A5,C5] = ntd(X, [10, 5, 90], opts);

% reconstruction
X_hat1 = ttm(C1, {A1{1}, A1{2}, A1{3}}, [1, 2, 3]);
X_hat3 = ttm(C3, {A3{1}, A3{2}, A3{3}}, [1, 2, 3]);
X_hat5 = ttm(C5, {A5{1}, A5{2}, A5{3}}, [1, 2, 3]);
%% Plots (Faces)
for i = 1:10
    figure;
    count = 1;
    % original
    for j = 1:9
        subplot(4, 9, count)
        imshow(reshape(double(X(i, j, :)), 192, 168), [])
        if j == 1; ylabel("original", "Interpreter","latex"); end
        count = count + 1;
    end
    % ilumination = 1
    for j = 1:9
        subplot(4, 9, count)
        imshow(reshape(double(X_hat1(i, j, :)), 192, 168), [])
        if j == 1; ylabel("ilumination = 1", "Interpreter","latex"); end
        count = count + 1;
    end
    % ilumination = 3
    for j = 1:9
        subplot(4, 9, count)
        imshow(reshape(double(X_hat3(i, j, :)), 192, 168), [])
        if j == 1; ylabel("ilumination = 3", "Interpreter","latex"); end
        count = count + 1;
    end
    % ilumination = 5
    for j = 1:9
        subplot(4, 9, count)
        imshow(reshape(double(X_hat5(i, j, :)), 192, 168), [])
        if j == 1; ylabel("ilumination = 5", "Interpreter","latex"); end
        count = count + 1;
    end
    sgtitle("Person " + num2str(i), "interpreter", "latex")
end
%% errors
errors = zeros(10, 9, 3);
for i = 1:10
    for j = 1:9
        errors(i, j, 1) = norm(X(i, j, :)-X_hat1(i, j, :));
        errors(i, j, 2) = norm(X(i, j, :)-X_hat3(i, j, :));
        errors(i, j, 3) = norm(X(i, j, :)-X_hat5(i, j, :));
    end
end
%% Plots (Errors)
figure;
for i = 1:10
    plot(errors(i, :, 1))
    hold on
end
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "interpreter", "latex", "NumColumns", 2)
title("Error for Each subject through different iluminations (components = 1)", "Interpreter","latex")
xlabel("ilumination", "Interpreter","latex")
ylabel("error", "Interpreter","latex")
grid minor

figure;
for i = 1:10
    plot(errors(i, :, 2))
    hold on
end
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "interpreter", "latex", "NumColumns", 2)
title("Error for Each subject through different iluminations (components = 3)", "Interpreter","latex")
xlabel("ilumination", "Interpreter","latex")
ylabel("error", "Interpreter","latex")
grid minor

figure;
for i = 1:10
    plot(errors(i, :, 3))
    hold on
end
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "interpreter", "latex", "NumColumns", 2)
title("Error for Each subject through different iluminations (components = 5)", "Interpreter","latex")
xlabel("ilumination", "Interpreter","latex")
ylabel("error", "Interpreter","latex")
grid minor

figure;
plot(mean(errors(:, :, 1), 2))
hold on
plot(mean(errors(:, :, 2), 2))
plot(mean(errors(:, :, 3), 2))
xlabel("subject", "Interpreter","latex")
ylabel("error")
title("Mean Error for each subject", "Interpreter","latex")
legend("components = 1", "components = 3", "components = 5", "interpreter", "latex")
grid minor
%% svd (part b)
T = zeros(32256, 90);
count = 1;
for i = 1:10
    for j = 1:9
        T(:, count) = double(X(i, j, :));
        count = count + 1;
    end
end
[U, S, V] = svd(T, "econ");
X_hat1_svd = U(:, 1:10) * S(1:10, 1:10) * V(:, 1:10)';
X_hat3_svd = U(:, 1:30) * S(1:30, 1:30) * V(:, 1:30)';
X_hat5_svd = U(:, 1:50) * S(1:50, 1:50) * V(:, 1:50)';
%% Plots (Faces)
countFace = 1;
for i = 1:10
    figure;
    count = 1;
    % original
    for j = 1:9
        subplot(4, 9, count)
        imshow(reshape(double(X(i, j, :)), 192, 168), [])
        if j == 1; ylabel("original", "Interpreter","latex"); end
        count = count + 1;
    end

   % 10 components
   for j = 1:9
       subplot(4, 9, count)
       imshow(reshape(X_hat1_svd(:, countFace), 192, 168), [])
       if j == 1; ylabel("components = 10", "Interpreter","latex"); end
       count = count + 1;
       countFace = countFace + 1;
   end

   % 30 components
   countFace = countFace - 9;
   for j = 1:9
       subplot(4, 9, count)
       imshow(reshape(X_hat3_svd(:, countFace), 192, 168), [])
       if j == 1; ylabel("components = 30", "Interpreter","latex"); end
       count = count + 1;
       countFace = countFace + 1;
   end

   % 50 components
   countFace = countFace - 9;
   for j = 1:9
       subplot(4, 9, count)
       imshow(reshape(X_hat5_svd(:, countFace), 192, 168), [])
       if j == 1; ylabel("components = 50", "Interpreter","latex"); end
       count = count + 1;
       countFace = countFace + 1;
   end
   sgtitle("Person " + num2str(i), "interpreter", "latex")
end
%% errors
X_hat1_svd_tensor = zeros(10, 9, 32256);
X_hat3_svd_tensor = zeros(10, 9, 32256);
X_hat5_svd_tensor = zeros(10, 9, 32256);
count = 1;
for i = 1:10
    for j = 1:9
        X_hat1_svd_tensor(i, j, :) = X_hat1_svd(:, count);
        X_hat3_svd_tensor(i, j, :) = X_hat3_svd(:, count);
        X_hat5_svd_tensor(i, j, :) = X_hat5_svd(:, count);
        count = count + 1;
    end
end

X_hat1_svd_tensor = tensor(X_hat1_svd_tensor);
X_hat3_svd_tensor = tensor(X_hat3_svd_tensor);
X_hat5_svd_tensor = tensor(X_hat5_svd_tensor);

errors_svd = zeros(10, 9, 3);
for i = 1:10
    for j = 1:9
        errors_svd(i, j, 1) = norm(X(i, j, :)-X_hat1_svd_tensor(i, j, :));
        errors_svd(i, j, 2) = norm(X(i, j, :)-X_hat3_svd_tensor(i, j, :));
        errors_svd(i, j, 3) = norm(X(i, j, :)-X_hat5_svd_tensor(i, j, :));
    end
end
%% Plots (Errors)
figure;
for i = 1:10
    plot(errors_svd(i, :, 1))
    hold on
end
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "interpreter", "latex", "NumColumns", 2)
title("Error for Each subject through different iluminations (components = 1)", "Interpreter","latex")
xlabel("ilumination", "Interpreter","latex")
ylabel("error", "Interpreter","latex")
grid minor

figure;
for i = 1:10
    plot(errors_svd(i, :, 2))
    hold on
end
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "interpreter", "latex", "NumColumns", 2)
title("Error for Each subject through different iluminations (components = 3)", "Interpreter","latex")
xlabel("ilumination", "Interpreter","latex")
ylabel("error", "Interpreter","latex")
grid minor

figure;
for i = 1:10
    plot(errors_svd(i, :, 3))
    hold on
end
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "interpreter", "latex", "NumColumns", 2)
title("Error for Each subject through different iluminations (components = 5)", "Interpreter","latex")
xlabel("ilumination", "Interpreter","latex")
ylabel("error", "Interpreter","latex")
grid minor

figure;
plot(mean(errors_svd(:, :, 1), 2))
hold on
plot(mean(errors_svd(:, :, 2), 2))
plot(mean(errors_svd(:, :, 3), 2))
xlabel("subject", "Interpreter","latex")
ylabel("error")
title("Mean Error for each subject", "Interpreter","latex")
legend("components = 1", "components = 3", "components = 5", "interpreter", "latex")
grid minor
%% Error Tucker vs. SVD
figure;
plot(mean(errors(:, :, 1), 2), "Color","blue")
hold on
plot(mean(errors_svd(:, :, 1), 2), "Color","red")
plot(mean(errors(:, :, 2), 2), "Color","blue")
plot(mean(errors(:, :, 3), 2), "Color","blue")
plot(mean(errors_svd(:, :, 2), 2), "Color","red")
plot(mean(errors_svd(:, :, 3), 2), "Color","red")
xlabel("subject", "Interpreter","latex")
ylabel("error")
title("Mean Error for each subject (Tucker vs. SVD)", "Interpreter","latex")
legend("Tucker", "SVD", "interpreter", "latex")
grid minor













