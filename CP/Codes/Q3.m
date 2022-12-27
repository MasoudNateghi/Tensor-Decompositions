clear; close all; clc;
load("amino.mat")
I = 5; 
J = 201;
K = 61;
T = zeros(I, J, K);
for k = 1:K
    T(:, :, k) = X(:, 1+(k-1)*J:k*J);
end

R = [2, 3, 4, 5];
for i = 1:length(R)
    U1_0 = randn(I, R(i));
    U2_0 = randn(J, R(i));
    U3_0 = randn(K, R(i));
    U0 = {U1_0, U2_0, U3_0};
    
    U = cpd_als(T,U0);
    U2 = U{2};
    figure(i);
    subplot(3, 2, 1)
    for j = 1:R(i)
        plot(250:450, U2(:, j))
        hold on
    end
    title("ALS", "Interpreter","latex")
    xlabel("Emmision Wavelength (nm)", "Interpreter","latex")
    ylabel("intensity", "Interpreter","latex")
    figure(length(R)+i)
    subplot(3, 2, 1)
    Consistency = corcond(T,U,[],1);
    title("Core consistency (ALS) "+ num2str(Consistency), "Interpreter","latex")

    [U1, U2, U3] = cp_als(T,U1_0, U2_0, U3_0);
    U = {U1, U2, U3};
    figure(i);
    subplot(3, 2, 2)
    for j = 1:R(i)
        plot(250:450, U2(:, j))
        hold on
    end
    title("ALS (developed)", "Interpreter","latex")
    xlabel("Emmision Wavelength (nm)", "Interpreter","latex")
    ylabel("intensity", "Interpreter","latex")
    figure(length(R)+i)
    subplot(3, 2, 2)
    Consistency = corcond(T,U,[],1);
    title("Core consistency (ALS (developed)) "+ num2str(Consistency), "Interpreter","latex")

    U = cpd3_sd(T,U0);
    U2 = U{2};
    figure(i);
    subplot(3, 2, 3)
    for j = 1:R(i)
        plot(250:450, U2(:, j))
        hold on
    end
    title("SD", "Interpreter","latex")
    xlabel("Emmision Wavelength (nm)", "Interpreter","latex")
    ylabel("intensity", "Interpreter","latex")
    figure(length(R)+i)
    subplot(3, 2, 3)
    Consistency = corcond(T,U,[],1);
    title("Core consistency (SD) "+ num2str(Consistency), "Interpreter","latex")

    U = cpd_gevd(T,R(i));
    U2 = U{2};
    figure(i);
    subplot(3, 2, 4)
    for j = 1:R(i)
        plot(250:450, U2(:, j))
        hold on
    end
    title("GEVD", "Interpreter","latex")
    xlabel("Emmision Wavelength (nm)", "Interpreter","latex")
    ylabel("intensity", "Interpreter","latex")
    figure(length(R)+i)
    subplot(3, 2, 4)
    Consistency = corcond(T,U,[],1);
    title("Core consistency (GEVD) "+ num2str(Consistency), "Interpreter","latex")

    U = cpd_minf(T,U0);
    U2 = U{2};
    figure(i);
    subplot(3, 2, 5)
    for j = 1:R(i)
        plot(250:450, U2(:, j))
        hold on
    end
    title("MINF", "Interpreter","latex")
    xlabel("Emmision Wavelength (nm)", "Interpreter","latex")
    ylabel("intensity", "Interpreter","latex")
    figure(length(R)+i)
    subplot(3, 2, 5)
    Consistency = corcond(T,U,[],1);
    title("Core consistency (MINF) "+ num2str(Consistency), "Interpreter","latex")

    U = cpd3_sgsd(T,U0);
    U2 = U{2};
    figure(i);
    subplot(3, 2, 6)
    for j = 1:R(i)
        plot(250:450, U2(:, j))
        hold on
    end
    title("SGSD", "Interpreter","latex")
    xlabel("Emmision Wavelength (nm)", "Interpreter","latex")
    ylabel("intensity", "Interpreter","latex")
    figure(length(R)+i)
    subplot(3, 2, 6)
    Consistency = corcond(T,U,[],1);
    title("Core consistency (SGSD) "+ num2str(Consistency), "Interpreter","latex")

    figure(i);
    sgtitle("Canonical Polyadic Tensor Decomposition (R = " + num2str(R(i)) + ")", "interpreter", "latex")

    figure(length(R)+i);
    sgtitle("Corcondia Criterion (R = " + num2str(R(i)) + ")", "interpreter", "latex")
end
