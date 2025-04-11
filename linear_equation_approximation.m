% Jacobi and Gauss-Seidel
% not yet
% Know that if A is diagonally dominant matrix, then both methods converge to the solution Ax = b.
% First, LU factorization..
clear;
A = [4 -1 -1 0;-1 4 0 -1;-1 0 4 -1;0 -1 -1 4];
b = [-1;3;7;11];
%LU factorization
function [L,U] = mylu(n,A)
%no need to initialize L and U since it is dynamically sized
    for k = 1:n-1
        for i = k+1:n
            A(i,k) = A(i,k)/A(k,k);
            for j=k+1:n
                A(i,j) = A(i,j)-A(i,k)*A(k,j);
            end
        end
    end
    L = tril(A,-1) + eye(n);
    U = triu(A);
    disp('L = ');
    disp(L);
    disp('U = ');
    disp(U);
end
function y = lowtri(n,L,b)
    y = zeros(n,1); %But, solving linear systems iteratively, MATLAB won't know the size of the vectors(from chatGPT)
    for k=1:n %forward
        y(k) = (b(k)-L(k,1:k-1)*y(1:k-1))/L(k,k);
    end
end

function x = uppertri(n,U,y)
    x = zeros(n,1);
    for k = n:-1:1 %backward
        x(k) = (y(k)-U(k,k+1:n)*x(k+1:n))/U(k,k);
    end
end
[L,U] = mylu(4,A);
y = lowtri(4,L,b);
x = uppertri(4,U,y);
disp('x = ');
disp(x);

