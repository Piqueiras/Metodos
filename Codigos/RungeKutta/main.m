clear all
fid = 1;

A = [(3+sqrt(3))/6, 0; -sqrt(3)/3, (3+sqrt(3))/6];
b = [1/2, 1/2];
c = [(3+sqrt(3))/6, (3-sqrt(3))/6];

% Definir la función f
f = @(x,Y) [Y(2);Y(2)*(Y(2)-1)/Y(1)];

% Solución exacta
sol = @(x) [(1+3*exp(-8*x))/8; 
            -3*exp(-8*x)];

% Parámetros del problema
a = 0;      %x0
d = 1;      %xN
eta = [1/2; -3];    %y0

% Tabla de errores
for h = [0.1, 0.01, 0.001]

    N = ceil((d-a) / h);

    % Resolver numéricamente con el método de Euler
    [xN, yN] = DIRK_yn(fid,f,a,eta,h,N,A,b,c,sol);
    
    % Mostrar resultados
    fprintf('h = %f, |y(xn) - yn| = %.10f, |y''(xn) - zn| = %.10f\n', h, abs(sol(xN(end))-yN(:,end)));
end

