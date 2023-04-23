function test6
% Autor: Miłosz Woźny, 320751
fprintf("Test sprawdza poprawność znajdowania przez funkcję\nzespolonych " + ...
    "pierwiastków macierzy trójdiagonalnych. Parametr mu podawany\ndo funkcji" + ...
    " jest zespolony w każdym z przykładów.\nObliczenia są wykonywane z domyślnymi" + ...
    " parametrami tol i max_it.\n\n")
pause;

A1 = [1, -1; 1, 1];
w1 = [1 + 1i, 1 - 1i];
mu1 = [0.5i, 2 - 3i];

A2 = [1, -1, 0, 0; 1, 1, -1, 0; 0, 1, 1, -1; 0, 0, 1, 1];
w2 = [1 + (1+sqrt(5))/2 * 1i, 1 - (1+sqrt(5))/2 * 1i, 1 + (sqrt(5) - 1)/2 * 1i, 1 - (sqrt(5) - 1)/2 * 1i];
mu2 = [0.5 + 3.5i, 2 - 2i, 1 + 0.25i, 0.75 - 0.99i];

A3 = [2 1 0 0 0; -1 2 1 0 0; 0 -1 2 1 0; 0 0 -1 2 1; 0 0 0 -1 2];
w3 = [2 + sqrt(3)*1i, 2 - sqrt(3)*1i, 2, 2 + 1i, 2-1i];
mu3 = [1 + 1.5i, 3 - 4i, 1 + 0.1i ,2.7 + 0.6i, 2.2 - sqrt(2)/2*1i];

fprintf("Macierz:\n");
disp(A1);
fprintf("|     mu     | w dokładna | w obliczona | błąd bezwzględny |\n");
for i = 1:2
    w = P2Z45_MWO_inverse_power_Givens(diag(A1, -1), diag(A1), diag(A1, 1), mu1(i));
    
    %fprintf("| %s | %s | %s | %16e |\n", real(mu1(i)), ...
    %    imag(mu1(i)), real(w1(i)), imag(w1(i)), real(w), imag(w), abs(w1(i) - w));
    fprintf("| %10s | %10s | %11s | %16e |\n", num2str(mu1(i)), num2str(w1(i)), num2str(w), abs(w1(i) - w));
end
pause;
fprintf("\n\nMacierz:\n");
disp(A2);
fprintf("|     mu     | w dokładna | w obliczona | błąd bezwzględny |\n");
for i = 1:4
    w = P2Z45_MWO_inverse_power_Givens(diag(A2, -1), diag(A2), diag(A2, 1), mu2(i));
    
    %fprintf("| %s | %s | %s | %16e |\n", real(mu1(i)), ...
    %    imag(mu1(i)), real(w1(i)), imag(w1(i)), real(w), imag(w), abs(w1(i) - w));
    fprintf("| %10s | %10s | %11s | %16e |\n", num2str(mu2(i)), num2str(w2(i)), num2str(w), abs(w2(i) - w));
end
pause;

fprintf("\n\nMacierz:\n");
disp(A3);
fprintf("|      mu      | w dokładna |  w obliczona  | błąd bezwzględny |\n");
for i = 1:5
    w = P2Z45_MWO_inverse_power_Givens(diag(A3, -1), diag(A3), diag(A3, 1), mu3(i));
    fprintf("| %12s | %10s | %13s | %16e |\n", num2str(complex(mu3(i))), num2str(complex(w3(i))), num2str(complex(w)), abs(w3(i) - w));
end



end

