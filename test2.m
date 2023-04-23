function test2
% Autor: Miłosz Woźny, 320751

fprintf("Funkcja sprawdza, jak zmienia się dokładność\nprzybliżenia wartości" + ...
    " własnej macierzy wraz z kolejnymi iteracjami\nmetody. Funkcja rysuje wykres" + ...
    " zależności błędu bezwzględnego\ndo liczby iteracji w skali logarytmicznej " + ...
    "oraz przybliża\nwspółczynnik zbieżności metody dla konkretnych danych.\n" + ...
    "Wartości własne są znane, a błąd wektora\nwłasnego otrzymywany jest przez" + ...
    "wcześniejsze obliczenie\ntego wektora z dużą dokładnością.\n\n")


A = zeros(4, 4, 8);
A(:, :, 1) = diag([1, 1.5, 3, 7]);
A(:, :, 2) = diag([1, 1.9, 4.7, 11]);
A(:, :, 3) = diag([1, 1.5, 3, 7]) + diag([2, 3.4, 5], 1);
A(:, :, 4) = diag([1, 5.6, 9.9, 13.49]) + diag([11, -7.4, 5], 1);
A(:, :, 5) = diag([2.24, 2.25, 2.4, 2.7]) + diag([1.4, 2, 6], -1);
A(:, :, 6) = diag([3, 5, 7, 8]) + diag([5.4, 2.1, -6.4], -1);
A(:, :, 7) = [1, 3, 0, 0; 2, 1, 5, 0; 0, 5, 1, 2; 0, 0, 3, 1];
A(:, :, 8) = [1, -2, 0, 0; -2, 1, 4, 0; 0, 4, -1, 2; 0, 0, 2, -1];

% i-ty wiersz to wartości własne i-tej macierzy z A
eig_values = [1, 1.5, 3, 7; 1, 1.9, 4.7, 11; 1, 1.5, 3, 7; 1, 5.6, 9.9, 13.49; ...
    2.24, 2.25, 2.4, 2.7; 7, 8, 3, 5; 0, 2, -5, 7; 5, 1, -1, -5];


mu = [0.1, 0.2, 0, 1.1, 0, 7.4, -0.6, 4.1];

for i = 1:8
    pause;
    [~,~,b_precise] = P2Z45_MWO_inverse_power_Givens(diag(A(:, :, i), -1), ...
        diag(A(:, :, i)), diag(A(:, :, i), 1), mu(i), 1e-16, 1000);
    [w, it, b] = P2Z45_MWO_inverse_power_Givens(diag(A(:, :, i), -1), ...
        diag(A(:, :, i)), diag(A(:, :, i), 1), mu(i), 1e-11, 500, 1);
    w_usable = w(1:it+1);
    b_usable = b(:, 1:it+1);
    err_w = abs(w_usable - eig_values(i, 1));
    err_w = err_w(err_w ~= 0);
    err_b = vecnorm(abs(b_usable) - abs(b_precise), 1);
    err_b = err_b(err_b ~= 0)';

    figure(2*i+1);
    semilogy(1:length(err_w), err_w);
    title("Błąd bezwzględny przybliżenia wartości własnej w zależności od ilości iteracji")
    xlabel("liczba iteracji");
    ylabel("logarytm z błędu bezwzględnego");

    figure(2*i+2);
    semilogy(1:length(err_b), err_b);
    title("Błąd bezwzględny przybliżenia wektora własnego w zależności od ilości iteracji")
    xlabel("liczba iteracji");
    ylabel("logarytm z błędu bezwzględnego");

    coeff_w = [ones(length(err_w) - 4, 1), (5:length(err_w))'] \ log10(err_w(5:end));
    fprintf("Szukana wartość własna: %f\n", eig_values(i, 1));
    fprintf("Parametr mu: %f\n", mu(i));
    fprintf("Druga najbliższa wartość własna: %f\n", eig_values(i, 2));
    fprintf("Dokładność rozwiązania: %e\n", err_w(end));
    fprintf("Ilość wykorzystanych iteracji: %d\n", it);
    fprintf("Oszacowany czynnik zbieżności wartości własnej w: %f\n", 10 ^ coeff_w(2));

    coeff_b = [ones(length(err_b) - 4, 1), (5:length(err_b))'] \ log10(err_b(5:end));
    fprintf("Oszacowany czynnik zbieżności wektora własnego b: %f\n\n", 10 ^ coeff_b(2));
end
end

