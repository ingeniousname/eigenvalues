function [d2_r, d3_r, d4_r, c, s, A] = Givens_rotate_matrix(d1, d2, d3, convert)
% Projekt 2, zadanie 45
% Miłosz Woźny, 320751
%
% Funkcja stosuje n-1 rotacji Givensa dla macierzy trójdiagonalnej (zadanej
% przez wektory wartości na przekątnych d1, d2, d3) w taki sposób, że po
% przekształceniu otrzymujemy macierz trójkątną górną, której wartości
% elementów zapisane są w wektorach d2_r, d3_r, d4_r (korzystając z 
% własności przekształcenia zauważam, że wszystkie niezerowe elementy tej
% macierzy znajdują się na przekątnych: głównej i dwóch kolejnych nad nią).
% Kolejne wartości cosinusów i sinusów kąta obrotu zastosowanego w
% kolejnych obrotach Givensa zapisywane są w wektorach c i s. Działanie
% funkcji jest analogiczne do znalezienia rozkładu D = QR macierzy, gdzie R
% jest macierzą trójkątną górną, a Q - macierzą unitarną, ale dane o
% przekształceniach zamiast w macierzy zapisane są w wektorach c i s.
% Funkcja zawiera dodatkowy parametr convert, kiedy jest on ustawiony na 1
% funkcja dodatkowo zwraca macierz A - reprezentacja macierzowa R 
% (domyślnie 0, niezalecane dla dużych rozmiarów macierzy).
% Wejście:
%       d1 - wektor elementów z dolnej przekątnej rozważanej macierzy 
%            trójdiagonalnej, długości n-1
%       d2 - wektor elementów ze środkowej (głównej) przekątnej rozważanej 
%            macierzy trójdiagonalnej, długości n
%       d3 - wektor elementów z górnej przekątnej rozważanej macierzy 
%            trójdiagonalnej, długości n-1
%       convert - opcjonalny parametr, ustawiony na 1 sprawia, że funkcja
%            dodatkowo zwraca R w postaci macierzy, domyślnie 0
% Wyjście:
%       d2_r - wektor elementów z głównej przekątnej zwracanej macierzy R,
%              długości n
%       d3_r - wektor elementów z górnej przekątnej zwracanej macierzy R,
%              długości n-1
%       d4_r - wektor elementów z przekątnej nad górną przekątną zwracanej 
%              macierzy R, długości n-2
%          c - wektor wartości cosinusów kątów kolejnych rotacji Givensa
%              zastosowanych na wierszach macierzy wejściowej
%          s - wektor wartości sinusów kątów kolejnych rotacji Givensa
%              zastosowanych na wierszach macierzy wejściowej
%          A - opcjonalnie, znaleziona macierz R zapisana w postaci 
%              macierzy
if nargin < 4
    convert = 0;
end

% inicjalizuje wektory
c = zeros(1,length(d1));
s = zeros(1,length(d1));
d2_r = d2;
d3_r = d3;
d4_r = zeros(length(d3)-1,1);

for i = 1:length(d1)
    % najpierw wykonuję rotację na płaszczyźnie odpowiadającej wierszom i
    % oraz i + 1 - znajdują się na nich elementy na głównej oraz dolnej 
    % przekątnej, wartość na dolnej przekątnej jest zerowana, a na górnej 
    % odpowiednio zmieniana
    

    r = hypot(d2_r(i), d1(i));
    if r > 0
        % zapisuję wartości sinusa i cosinusa kąta, o jaki trzeba obrócić
        % wektor, aby otrzymać w macierzy zero na dolnej przekątnej
        c(i) = d2_r(i)'/r;
        s(i) = -d1(i)'/r;
        d2_r(i) = r;

        % następnie wykonuję analogiczną rotację dla pozostałych wektorów
        % kolumnowych w macierzy, które mają niezerowe elementy w wierszach
        % i-tym oraz (i+1)-szym, jest to macierz trójdiagonalna więc istnieją
        % tylko dwa takie wektory
        u1 = Givens_rotate_vector2([d3_r(i) d2_r(i + 1)],c(i),s(i));
        d3_r(i) = u1(1);
        d2_r(i+1) = u1(2);

        if i < length(d1)
            u2 = Givens_rotate_vector2([0 d3(i + 1)],c(i),s(i));
            d4_r(i) = u2(1);
            d3_r(i+1) = u2(2);
        end
    end
end


% opcjonalna konwersja wynikowej macierzy R do macierzy
if convert == 1
    A = eye(length(d2)).*d2_r;
    A(1, 2) = d3_r(1);
    for i = 1:length(d4_r)
        A(i+1,i+2) = d3_r(i+1);
        A(i,i+2) = d4_r(i);
    end
end

end % function

