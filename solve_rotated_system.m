function x = solve_rotated_system(d2, d3, d4, y)
% Projekt 2, zadanie 45
% Miłosz Woźny, 320751
%
% Funkcja rozwiązuje układ równań w postaci macierzy trójkątnej górnej 
% zadanej przez trzy diagonale d2, d3, d4 (macierzy zwróconej przez funkcję 
% Givens_rotate_matrix)
% Wejście:
%       d2_r - wektor elementów z głównej przekątnej zwracanej macierzy R,
%              długości n
%       d3_r - wektor elementów z górnej przekątnej zwracanej macierzy R,
%              długości n-1
%       d4_r - wektor elementów z przekątnej nad górną przekątną zwracanej 
%              macierzy R, długości n-2
%          y - wektor z równania Rx = y
% Wyjście:
%          x - wektor długości length(y), rozwiązanie układu równań Rx = y

n = length(d2);
% inicjuję wektor rozwiązań
x = zeros(n,1);
% bezpośrednio obliczam dwa pierwsze rozwiązania
x(n) = y(n)/d2(n);
x(n-1) = (y(n-1)-d3(n-1)*x(n))/d2(n-1);

% iteracyjnie obliczam pozostałe rozwiązania
for i = n-2:-1:1
    x(i) = (y(i)-d3(i)*x(i+1)-d4(i)*x(i+2))/d2(i);
end

end % function