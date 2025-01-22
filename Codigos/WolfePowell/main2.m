% Definimos la función de costo J(x) y su gradiente
J = @(x) (x(1)-1)^2 + (x(2)-x(1)^3)^2;
dJdx = @(x) [2*(x(1)-1) - 6*(x(2)-x(1)^3)*x(1)^2;
             2*(x(2)-x(1)^3)];

% Parámetros iniciales
x0 = [1; 2];  % Punto inicial
tol = 1e-6;   % Tolerancia
itmax = 1000; % Número máximo de iteraciones
rho = 0.2;    % Parámetro Wolfe
sigma = 0.4;  % Parámetro Wolfe
fid = 1;      % Salida a consola (1 para pantalla)

% Ejecutar el método de gradiente con Wolfe
[x_opt, nor_grad, iter] = GradWolfeConj(J, dJdx, x0, itmax, tol, rho, sigma, fid);

% Mostrar resultados
fprintf('\nResultados:\n');
fprintf('Iteraciones: %d\n', iter);
fprintf('Norma del gradiente final: %.6e\n', nor_grad);
fprintf('Solución aproximada: x_1 = %.6f, x_2 = %.6f\n', x_opt(1), x_opt(2));

% Crear tabla de resultados
fprintf('\n|  k  |  ||g_k||  |  (x_k)_1  |  (x_k)_2  |\n');
fprintf('|----|-----------|-----------|-----------|\n');
fprintf('| %2d | %.6e | %.6f | %.6f |\n', iter, nor_grad, x_opt(1), x_opt(2));