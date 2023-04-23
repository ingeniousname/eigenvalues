function [w, it, b] = P2Z45_MWO_inverse_power_Givens(d1, d2, d3, mu, tol, max_it, series)
% Projekt 2, zadanie 45
% Miłosz Woźny, 320751
%
% Funkcja szacuje wartość własną macierzy trójdiagonalnej D (o elementach 
% na przekątnych zdefiniowanych w wektorach d1, d2, d3) najbliższą
% wartości mu przy użyciu odwrotnej metody potęgowej z normowaniem
% zastosowanej dla macierzy D - mu*I, której rozkład, wykorzystywany przy 
% rozwiązywaniu układów równań, znajdowany jest przy pomocy 
% obrotów Givensa. Jeżeli po przekroczeniu maksymalnej ilości
% iteracji max_it określona dokładność tol nie zostanie osiągnięta,
% funkcja zwraca najlepsze obecne przybliżenie, w przeciwnym wypadku 
% zwraca przybliżoną wartość szukanej wartości własnej macierzy D.
% Funkcja szuka zespolonych wartości własnych wtedy i tylko wtedy, gdy
% parametr mu jest liczbą zespoloną.
% Wejście:
%       d1 - wektor pionowy elementów z dolnej przekątnej rozważanej 
%            macierzy trójdiagonalnej, długości n-1, o elementach 
%            rzeczywistych
%       d2 - wektor pionowy elementów ze środkowej (głównej) przekątnej 
%            rozważanej macierzy trójdiagonalnej, długości n, o elementach 
%            rzeczywistych
%       d3 - wektor pionowy elementów z górnej przekątnej rozważanej 
%            macierzy trójdiagonalnej, długości n-1, o elementach 
%            rzeczywistych
%       mu - parametr określający, którą wartość własną macierzy D
%            próbujemy znaleźć (najbliższą wartości mu), musi być on różny
%            niż któraś z wartości własnych macierzy, może być zespolony
%       tol - parametr określający dokładność szukanego rozwiązania,
%            występuje w warunku stopu, domyślnie 1e-7
%       max_it - liczba naturalna, maksymalna ilość iteracji, domyślnie 100
%       series - opcjonalny parametr, ustawiony na 1 sprawia, że w nie jest
%       pojedynczą wartością, a wektorem kolejnych przybliżeń, a b nie jest
%       ostatnim przybliżeniem wektora własnego, a macierzą, której i-ta
%       kolumna jest i-tym obliczonym przybliżeniem wektora własnego
% Wyjście:
%       w - wartość liczbowa szukanej wartości własnej, lub w przypadku
%           niepowodzenia (przekroczenia ilości iteracji) - obecne
%           przybliżenie
%       it - ilość przeprowadzonych iteracji
%       b - ostatnie przybliżenie wektora własnego macierzy

% ustawiam domyślne parametry
if nargin < 7
    if nargin < 6
        if nargin < 5
           tol = 1e-7; 
        end
        max_it = 100;
    end
    series = 0;
end

n = length(d2);
% Przypadek macierzy 1x1 (mało interesujący)
if n == 1
    w = d2;
    it = 0;
    b = 1;
    return;
end

% Znajduję rozkład macierzy D - mu*I = QR za pomocą rotacji Givensa, macierz zapisuję 
% w postaci trzech wektorów - diagonali otrzymanej macierzy oraz w
% wektorach c i s zapisuje wartości kolejno cosinusów i sinusów kątów 
% kolejnych obrotów (inna reprezentacja danych zawartych przez macierz Q)

[d2_r,d3_r,d4_r,c,s] = Givens_rotate_matrix(d1,d2-mu,d3);

% Inicjuję wektor własny losowymi wartościami i go normuję (inicjuję wektor
% jako zespolony wtedy i tylko wtedy gdy mu jest zespolone)
if isreal(mu)
    b = rand(n, 1);
else
    b = complex(rand(n, 1), rand(n, 1));
end
b = b/norm(b,2);
b_next = b;

% inicjowanie rozwiązania, jeżeli flaga series ustawiona jest na 1
if series == 1
    w = zeros(max_it + 1, 1);
    w(1) = b'*tridiagonal_product(b,d1,d2,d3);
    b_ret = zeros(n, max_it + 1);
    b_ret(:, 1) = b;
end

% Pętla iteracyjna
for i = 1:max_it
    b = b_next;
    % obliczam kolejne przybliżenie wektora własnego macierzy D, zgodnie z
    % równaniem b_{k+1} = (QR)^-1 * b_k = R^-1 * (Q'*b_k) = R \ (Q'*b_k),
    % najpierw stosuję na wektorze b_k odpowiednie przekształcenie, a
    % następnie rozwiązuję układ równań
    b_next = Givens_rotate_vector(b_next,c,s);
    b_next = solve_rotated_system(d2_r,d3_r,d4_r,b_next);
    % normuję wektor b_{k+1}
    b_next = b_next/norm(b_next,2);
    % jeżeli series jest ustawione na 1 - zapisuje wyniki z obecnej
    % iteracji do wektora/macierzy
    if series == 1
        w(i+1) = b_next'*tridiagonal_product(b_next,d1,d2,d3);
        b_ret(:, i+1) = b_next;
    end
    % sprawdzam warunek stopu
    idx_kp1 = get_max_idx(b_next);
    idx_k = get_max_idx(b);
    if norm(abs(b_next(idx_kp1))/b_next(idx_kp1)*b_next-...
            abs(b(idx_k))/b(idx_k)*b,2) < tol
        b = b_next;
        break;
    end
    
end

it = i;

% zwracam wartość własną, obliczam ją według wzoru b' * D * b, gdzie b jest
% wektorem własnym D, nie muszę dzielić przez b'*b ponieważ stosuję normę
% drugą przy normowaniu 
if series == 0
    w = b_next'*tridiagonal_product(b_next,d1,d2,d3);
else
    w(it+1) = b_next'*tridiagonal_product(b_next,d1,d2,d3);
    b = b_ret;
end

end % function

