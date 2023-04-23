function test5
% Autor: Miłosz Woźny, 320751

fprintf("Test sprawdza poprawność funkcji dla macierzy trójdiagonalnych\no " + ...
    "dużych rozmiarach. Sprawdzana jest wartość błędu:\n|Ax - wx|\ndla znalezionej" + ...
    " wartości własnej w najbliższej zeru.\nMacierze są losowane z ustalonym ziarnem, a obliczenia\ndokonywane są z " + ...
    "dokładnością tol = 1e-11, oraz z maksymalną liczbą\niteracji max_it = 500.\n\n")
rng(1);
pause;

n_vector = [1000 2000 10000 50000 100000 200000];
for n = n_vector
    d1 = 1000 * rand(n - 1, 1);
    d2 = 1000 * rand(n, 1);
    d3 = 1000 * rand(n - 1, 1);
    tic
    [w, ~, b] = P2Z45_MWO_inverse_power_Givens(d1, d2, d3, 0, 1e-11, 500);
    t = toc;
    err = norm(tridiagonal_product(b, d1, d2, d3) - w*b,1);
    fprintf("Błąd |Ax - wx| dla macierzy %dx%d: %e\n", n, n, err);
    fprintf("Czas trwania obliczeń: %fs\n\n", t);
end


end

