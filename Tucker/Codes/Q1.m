clear; close all; clc;
X = rand(4, 3, 2);
[G, U1, U2, U3] = hooi(X, 3, 2, 2);
X_hat = kruskal_tucker(G, U1, U2, U3);
tensor_norm(X_hat-X)