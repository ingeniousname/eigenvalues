 A = [1 3 0 0; 8 2 2 0; 0 9 5 1; 0 0 3 4];
 A_base = A;
 r1 = sqrt(A(1, 1)*A(1, 1) + A(2, 1)*A(2,1));
 Q1 = eye(4);
 Q1(1,1) = A(1,1) / r1;
 Q1(2,2) = A(1,1) / r1;
 Q1(2,1) = -A(2,1) / r1;
 Q1(1,2) = A(2,1) / r1;
    
 A = Q1 * A;
 r1 = sqrt(A(2, 2)*A(2, 2) + A(3, 2)*A(3,2));
 Q2 = eye(4);
 Q2(2,2) = A(2,2) / r1;
 Q2(3,3) = A(2,2) / r1;
 Q2(3,2) = -A(3,2) / r1;
 Q2(2,3) = A(3,2) / r1;
 
 A = Q2 * A;
 r1 = sqrt(A(3, 3)*A(3, 3) + A(4, 3)*A(4,3));
 Q3 = eye(4);
 Q3(3,3) = A(3,3) / r1;
 Q3(4,4) = A(3,3) / r1;
 Q3(4,3) = -A(4,3) / r1;
 Q3(3,4) = A(4,3) / r1;
    

%A = Q3 * A;
y = [1; 1; 1; 1];
[d2,d3,d4,c,s,R] = Givens_rotate_matrix(diag(A_base, -1), diag(A_base), diag(A_base, 1), 1);

%A_base*y
%Givens_rotate_vector_long(R*y, c, s, 'inverse')


%Q3*Q2*Q1*A_base*y

% Q*(A*y) = R*y
%R*y
%Givens_rotate_vector(A_base*y, c, s)

% A*y = Q'*(R*y)
%A_base*y
%Givens_rotate_vector(R*y, c, s, 'inverse')

A_base \ y
solve_rotated_system(d2, d3, d4, Givens_rotate_vector(y, c, s))
 


 


