function test3
% Autor: Miłosz Woźny, 320751
fprintf("Test sprawdza poprawność rozkładu A = QR otrzymanego za pomocą\nfunkcji " + ...
    "Givens_rotate_matrix. Obliczana jest norma wyrażenia |A - QR| dla\nstu losowych" + ...
    " macierzy trójdiagonalnych A rozmiaru 100x100. Rozważana norma\npowinna być bardzo mała.\n");

n = 100;
k = 100;


error = 0;
for it = 1:k
    d1 = rand(n-1, 1);
    d3 = rand(n-1, 1);
    d2 = rand(n, 1);
    [d2_r,d3_r,d4_r,c,s,R] = Givens_rotate_matrix(d1,d2,d3,1);
    A = diag(d1, -1) + diag(d2) + diag(d3, 1);
    QR = R;
    for i = 1:n
        QR(:, i) = Givens_rotate_vector(QR(:, i), c, s, 'inverse');
    end
    error = error + norm(A - QR, 'fro');

end

fprintf("Średnia norma wyrażenia |A - QR|: %e\n", error / k);

end


