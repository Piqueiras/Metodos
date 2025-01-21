clear all
fid = 1;

% y' = z
% z' = cos(3x) - y
% Definir la función f
f = @(x,Y) [Y(2);cos(3*x)-Y(1)];

% Solución exacta
sol = @(x) [(9*cos(x) - cos(3*x))/8; 
            (-9*sin(x) + 3*sin(3*x))/8];

% Parámetros del problema
a = 0;      %x0
b = 2.4;      %xN
eta = [1; 0];    %y0

% Tabla de errores
for h = [0.08, 0.04, 0.02, 0.01, 0.005]

    N = ceil((b-a) / h);

    % Resolver numéricamente con el método de Euler
    [xN, yN] = AdamsBashforth(fid,f,a,eta,h,N,sol);
    
    % Mostrar resultados
    fprintf('h = %f, |y(xn) - yn| = %.10f, |y''(xn) - zn| = %.10f\n', h, abs(sol(xN)-yN));
end

