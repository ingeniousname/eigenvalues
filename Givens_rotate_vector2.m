function u = Givens_rotate_vector2(v, c, s)
% Projekt 2, zadanie 45
% Miłosz Woźny, 320751
%
% Funkcja stosuję rotację Givensa dla wektora dwuwymiarowego, obracając go
% o kąt zadany wartością jego cosinusa c oraz sinusa s w kierunku
% przeciwnym do ruchu wskazówek zegara
% Wejście:
%       v - przekształcany wektor pionowy, length(v) = 2
%       c - wartość cosinusa kąta obrotu
%       s - wartość sinusa kąta obrotu
% Wyjście:
%       u - wektor v po przekształceniu

u(1) = c*v(1)-s*v(2);
u(2) = s'*v(1)+c'*v(2);

end % function

