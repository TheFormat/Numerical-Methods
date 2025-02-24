clear;
function [t, w] = myeuler(a, b, y0, N) % 오일러 방법
    h = (b - a) / N; % 간격
    t = zeros(1, N + 1);  % Array 모든 값이 0인 행렬
    w = zeros(1, N + 1);  % Array
    t(1) = a; %처음 element가 a
    w(1) = y0;
    % function 예제가 y' = y-t^2+1, 이 함수만 바꾸면 비선형 미분방정식도 작은 오차로 유추 가능
    f = @(t, w)w - t^2 + 1;
    for i = 1:N
        t(i + 1) = t(i) + h;
        w(i + 1) = w(i) + h * f(t(i), w(i));
    end
end

function [t,w] = myrk2(a,b,y0,N) % 2차 룽가쿠타 중에서도 modified Euler
    h = (b - a)/N;
    t = zeros(1,N + 1);
    w = zeros(1,N + 1);
    t(1) = a;
    w(1) = y0;
    f = @(t,w)w-t^2 + 1;
    for i = 1:N
        t(i + 1) = t(i) + h;
        w(i + 1) = w(i) + h*(f(t(i),w(i))/2 + f(t(i)+h,w(i)+h*f(t(i),w(i)))/2);
    end
end

%4차 Runge_Kutta는 추후에
%오차들도 업데이트 필요

a = 0; %0<t<2
b = 2;
y0 = 0.5; %IV가 0.5
N = 10; % 간격인 h는 0.2가 될 것

[teu,weu] = myeuler(a,b,y0,N); %함수를 만들었으니 넣어보자
[trk2, wrk2] = myrk2(a,b,y0,N);
rweu = round(weu,7); % array의 elements 모두 7자리 반올림
rwrk2 = round(wrk2,7);
y = @(t) (t+1).^2-0.5*exp(t); %exact solution
tval = 0:0.2:2; %array에 대입하기 위해서 t값들을 일일이 배정
y_exact = y(tval);
ry_exact = round(y_exact,7);
disp('Euler method');
disp('Solution values (w):');
disp(rweu);
disp('second order Runge-Kutta method')
disp('Solution values (w):');
disp(rwrk2);
disp('Exact value y(t):');
disp(ry_exact);
err_eu = abs(y_exact-weu); % error 오차
err_rk2 = abs(y_exact-wrk2);
disp('Absolute error of Euler Method')
disp(err_eu);
disp('Absolute error of Runge-Kutta')
disp(err_rk2);
%second order Runga-Kutta Method is more accurate than Euler method.
