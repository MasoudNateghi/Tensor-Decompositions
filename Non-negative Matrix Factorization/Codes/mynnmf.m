function [B, C] = mynnmf(A, j, algorithm, B0, C0, MaxItter, Tol)
n = size(A, 1);
m = size(A, 2);
if nargin == 2
    algorithm = "als";
    B0 = rand(n, j);
    C0 = rand(j, m);
elseif nargin == 5
    MaxItter = 100;
    Tol = 1e-4;
elseif nargin == 6
    Tol = 1e-4;
end
epsilon = 1e-16;
Itter = 1;
E = ones(size(A));
B = B0;
C = C0;
while Itter <= MaxItter && norm(E, "fro") >= Tol
    if algorithm == "als"
        C = max(epsilon, pinv(B' * B) * B' * A);
        B = max(epsilon, A * C' * pinv(C * C'));
    elseif algorithm == "mult"
        C = C .* (B' * A) ./ (B' * B * C + epsilon);
        B = B .* (A * C') ./ (B * C * C' + epsilon); %#ok<MHERM> 
    end
    E = A - B * C;
    Itter = Itter + 1;
end