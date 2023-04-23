function test4
% Autor: Miłosz Woźny, 320751

fprintf("Test sprawdza czas działania pojedynczej iteracji\nmetody w zależności" + ...
    " od rozmiaru danych.\nW teście zasymulowana jest pojedyncza iteracja z losowo\n" + ...
    "wygenerowanymi danymi. Wykres czasu od rozmiaru danych\npowinien być liniowy.")

tx = zeros(1000, 1);
n_vector = round(linspace(3, 10000, 1000));
for i = 1:1000
    n = n_vector(i);
    d2 = rand(n, 1);
    d3 = rand(n - 1, 1);
    d4 = rand(n - 2, 1);
    b = rand(n, 1);
    c = rand(n-1, 1);
    s = rand(n-1,1);
    tic;
    b_next = Givens_rotate_vector(b,c,s);
    b_next = solve_rotated_system(d2,d3,d4,b_next);
    b_next = b_next/norm(b,2);

    t = toc;
    tx(i) = t;
end

plot(1:1000, tx);
ylabel("czas, s");
xlabel("n");
title("Wykres zależności czasu działania iteracji od rozmiaru danych");

end

