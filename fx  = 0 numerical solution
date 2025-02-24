% Applied Numerical Methods(MATH-5336-A)
% DonghyoLee(Fred) EagleID: 901465856
% 11th, September, 2024
f = @(x)x^3+3*x^2-1; % 함수 예제
% Goal: find zeros of f(x) = 0. 
% Bisection Method
a0 = 0;
b0 = 1;
fa = f(a0);
fb = f(b0);
if (sign(fa) == sign(fb))
    fprintf('nope');% Just in case.. Intermediate Value Theorem 사용
end
eps = 10^(-6); % 입실론 설정
tol = 10^(-6); % tolerance 설정, 허용오차
itcount = 0;
e = b0 -a0; % 간격이 되겠다
M = 50; % 시행 횟수
while(itcount<=M && e >= tol)
    itcount = itcount +1;
    e = e/2;% Every time the loop goes on, interval shorten by half
    p = a0 + e;
    fp = f(p);
    if (abs(fp)<eps) % 만약에 입실론보다 함숫값이 작아진다면 끝, 물론 잘 일어나진 않더라
        fprintf('%g, %g, %g',itcount,p,fp);
        return;
    end
    if (sign(fa)==sign(fp)) % 부호가 같으면 a와 p를 비교
        a0 = p;
        fa = fp;
        fprintf('%g, %g, %g, %g \n',a0,b0,p,fp);
    else %부호가 다르다면 b와 p를 비교
        b0 = p;
        fb = fp;
        fprintf('%g, %g, %g, %g \n',a0,b0,p,fp);
    end
end
fprintf('Bisection Method, after %g steps, zero is %g \n', itcount,p);

% Newton's Method, Fixed-Point iteration의 가장 보편적인 예
% 만약에 왜 뉴턴 방법으로 해를 구할 수 있는지는 Theorem으로 남김
x0 = 1;
df = @(x)3*x^2+6*x; % f(x) = 0 과 x = g(x) = x-f/f'의 해가 같으니까
v = f(x0);
Nitcount = 0;
while(Nitcount<=M)
    Nitcount = Nitcount +1;
    q = x0 - v/df(x0);
    v = f(q);
    fprintf('%g, %g\n',Nitcount,q);
    if (abs(q-x0) < tol || abs(v) < eps) % ||은 or연산자
        fprintf('Newton Method, after %g steps, zero is %g',Nitcount,q);
        return;% End of loop if following condition is satisfied
    end
    x0 = q;
end
