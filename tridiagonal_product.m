function u = tridiagonal_product(v, d1, d2, d3)
% Projekt 2, zadanie 45
% Miłosz Woźny, 320751
%
% Funkcja oblicza iloczyn A * v zwraca go w zmiennej u, przy czym 
% macierz A ma niezerowe elementy na trzech przekątnych, zapisana jest 
% w postaci trzech wektorów (zgodnie z wartością zwracaną funkcji 
% Givens_rotate_matrix)
% Wejście:
%       v - wektor pionowy, rozmiaru n
%       d1, d2, d3 - kolejno dolna, główna i górna przekątna
%           macierzy trójdiagonalnej wymiaru n x n
% Wyjście:
%       u - obliczony iloczyn

n = length(d2);

u = v.*d2;
u(1:n-1) = u(1:n-1)+d3.*v(2:n);
u(2:n) = u(2:n)+d1.*v(1:n-1);


end % function

