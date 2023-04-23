clear;
clc;
A = [1, -1; 1, 1];
[V, D] = eig(A);
v = V(:, 1);
d1 = [4;3;2;1];
d2 = [5;6;7;-2;-1];
d3 = [-1;-1;-1;-1];

A = diag(d1, -1) + diag(d2) + diag(d3, 1);
eig(A)

[d2_r, d3_r, d4_r, c, s, R] = Givens_rotate_matrix(1, [1; 1], -1, 1);
[w, it, b] = P2Z45_MWO_inverse_power_Givens(d1, d2, d3, 5, 1e-15, 10000);
err = norm(b*w - tridiagonal_product(b, d1, d2, d3),2);
%x = complex(rand(n, 1), rand(n, 1));
x = v;
x = x / norm(x, 2);



