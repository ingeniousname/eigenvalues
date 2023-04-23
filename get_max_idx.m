function i = get_max_idx(v)
% Projekt 2, zadanie 45
% Miłosz Woźny, 320751
%
% Funkcja znajduje indeks największego co do modułu elementu wektora v i go
% zwraca.

i = 1;
% znajduję największy co do modułu element w wektorze v
for j = 2:length(v)
    if abs(v(i)) < abs(v(j))
        i = j;
    end
end

end % function

