 %A = [1 3 0 0; 8 2 2 0; 0 9 5 1; 0 0 3 4];
 A = [1, 3, 0, 0; 2, 1, 5, 0; 0, 5, 1, 2; 0, 0, 3, 1];
 A_base = A;

n = 5;

%[w, it, b] = P2Z45_MWO_inverse_power_Givens(zeros(n-1, 1),[4.4455;4.4456;4.4457;4.4458;4.4459],zeros(n-1,1),4, 1e-9, 4.4458);

[w, it, b] = P2Z45_MWO_inverse_power_Givens(diag(A, -1),diag(A), diag(A, 1), -0.6, 1e-9, 500);


%[w, it, b] = P2Z45_MWO_inverse_power_Givens([1], [1;1], [-1], 0, 1e-11, 500, 1);

%k = power(10, 8);
%d12 = (4.4 + (double(1:5) ./ k))';
%v = rand(n, 1);
%w = v' * tridiagonal_product(v, zeros(4, 1), d12, zeros(4, 1));
%[w, it, b] = P2Z45_MWO_inverse_power_Givens(zeros(4, 1), d12, zeros(4, 1), 0, 1e-7, 1000, 1);

