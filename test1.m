 function test1
% Autor: Miłosz Woźny, 320751
fprintf("Test sprawdza dokładność rozwiązania zagadnienia własnego\n" + ...
    "dla macierzy diagonalnych przy użyciu parametrów domyślnych.\n" + ...
    "W teście wykorzystywane są macierze " + ...
    "diagonalne,\nponieważ w przeciwnym wypadku bardzo ciężko\n" + ...
    "otrzymać dokładne rozwiązanie do porównania z wynikiem funkcji.\n" + ...
    "w - wartość własna\n\n")


d1 = zeros(999, 1);
d12 = (1:1000)';
d3 = zeros(999, 1);

mu11 = 0.1;
mu12 = pi;
mu13 = 674.278;
mu14 = 774.12;
mu15 = 999.5;
pause;

fprintf("Macierz 1, Wektor na przekątnej: (1:1000)'\n");
fprintf("|    mu    | w dokładna | w obliczona | błąd bezwzględny |\n");
w = P2Z45_MWO_inverse_power_Givens(d1,d12,d3, mu11);

fprintf("| %8.4f | %10.4f | %11.4f | %-16e |\n", mu11, round(mu11) + 1, w, abs(round(mu11) + 1 - w));
w = P2Z45_MWO_inverse_power_Givens(d1,d12,d3, mu12);


fprintf("| %8.4f | %10.4f | %11.4f | %-16e |\n", mu12, round(mu12), w, abs(round(mu12) - w));
w = P2Z45_MWO_inverse_power_Givens(d1,d12,d3, mu13);

fprintf("| %8.4f | %10.4f | %11.4f | %-16e |\n", mu13, round(mu13), w, abs(round(mu13) - w));
w = P2Z45_MWO_inverse_power_Givens(d1,d12,d3, mu14);

fprintf("| %8.4f | %10.4f | %11.4f | %-16e |\n", mu14, round(mu14), w, abs(round(mu14) - w));
w = P2Z45_MWO_inverse_power_Givens(d1,d12,d3, mu15);

fprintf("| %8.4f | %10.4f | %11.4f | %-16e |\n", mu15, round(mu15), w, abs(round(mu15) - w));
pause;

d21 = zeros(199, 1);
d22 = ((-100:99) .^ 3)';
d23 = zeros(199, 1);
mu21 = 0.22;
mu22 = -23.334;
mu23 = 512;
mu24 = 9830.12318;
mu25 = -90744;
w21 = 0;
w22 = -27;
w23 = 512;
w24 = 9261;
w25 = -91125;
fprintf("\nMacierz 2, Wektor na przekątnej: ((-500:499) .^ 3)'\n");
fprintf("|     mu     |  w dokładna  |  w obliczona  | błąd bezwzględny |\n");
w = P2Z45_MWO_inverse_power_Givens(d21,d22,d23, mu21);

fprintf("| %10.3f | %12.4f | %13.4f | %16e |\n", mu21, w21, w, abs(w21 - w));
w = P2Z45_MWO_inverse_power_Givens(d21,d22,d23, mu22);

fprintf("| %10.3f | %12.4f | %13.4f | %16e |\n", mu22, w22, w, abs(w22 - w));
w = P2Z45_MWO_inverse_power_Givens(d21,d22,d23, mu23);

 fprintf("| %10.3f | %12.4f | %13.4f | %16e |\n", mu23, w23, w, abs(w23 - w));
w = P2Z45_MWO_inverse_power_Givens(d21,d22,d23, mu24);

fprintf("| %10.3f | %12.4f | %13.4f | %16e |\n", mu24, w24, w, abs(w24 - w));
w = P2Z45_MWO_inverse_power_Givens(d21,d22,d23, mu25);

fprintf("| %10.3f | %12.4f | %13.4f | %16e |\n", mu25, w25, w, abs(w25 - w));
end % function