function u = Givens_rotate_vector(v, c, s, option)
% Projekt 2, zadanie 45
% Miłosz Woźny, 320751
%
% Funkcja przeprowadza odpowiednie rotacje Givensa (zastosowane wcześniej
% na macierzy, zapisane w postaci wektorów wartości sinusów i cosinusów 
% odpowiednio s i c) na wektorze v, opcjonalnie opcja daje możliwość 
% wykonania operacji odwrotnej. Wywołanie funkcji jest równoważne, przy
% znanym już rozkładzie macierzy D = QR, wykonaniu mnożenia Qv (Q'v).
% Wejście:
%       v - przekształcany wektor
%       c - wektor cosinusów kolejnych przekształceń macierzy D, otrzymany
%           przez wywołanie funkcji Givens_rotate_matrix
%       s - wektor sinusów kolejnych przekształceń macierzy D, otrzymany
%           przez wywołanie funkcji Givens_rotate_matrix
%       option - opcjonalny parametr, ustawiony na 'inverse' dokonuje
%           przekształcenia odwrotnego, domyślnie 'default'
% Wyjście:
%       u - wektor v po przekształceniu


if nargin < 4
    option = 'default';
end

iter = 1:length(v)-1;

% Jeżeli wykonywane jest przekształcenie odwrotne, odwracam wektor iter 
% oraz neguję wartości sinusów, przy przekształceniu odwrotnym (u = Q'v) 
% obrotów należy dokonać w odwrotnej kolejności i o przeciwny kąt
if strcmp(option, 'inverse')
    iter = flip(iter);
    s = -s;
end

u = v;
for i = iter
    % wykonuję rotację na kolejnych płaszczyznach
    u(i:i+1) = Givens_rotate_vector2(u(i:i+1),c(i),s(i));
end

end

