% Definimos la función de costo J(x) y su gradiente
J = @(x) 2*x(1)^2 - 1.05*x(1)^4 + (1/6)*x(1)^6 + x(1)*x(2) + x(2)^2;
dJdx = @(x) [4*x(1) - 4.2*x(1)^3 + x(1)^5 + x(2);
             x(1) + 2*x(2)];

% Parámetros iniciales
x0 = [1; 2];  % Punto inicial
tol = 1e-8;   % Tolerancia
itmax = 1000; % Número máximo de iteraciones
rho = 0.2;    % Parámetro Wolfe
sigma = 0.4;  % Parámetro Wolfe
fid = 1;      % Salida a consola (1 para pantalla)

% Ejecutar el método de gradiente con Wolfe
[x_opt, nor_grad, iter] = GradWolfe(J, dJdx, x0, itmax, tol, rho, sigma, fid);

% Mostrar resultados
fprintf('\nResultados:\n');
fprintf('Iteraciones: %d\n', iter);
fprintf('Norma del gradiente final: %.6e\n', nor_grad);
fprintf('Solución aproximada: x_1 = %.6f, x_2 = %.6f\n', x_opt(1), x_opt(2));

% Crear tabla de resultados
fprintf('\n|  k  |  ||g_k||  |  (x_k)_1  |  (x_k)_2  |\n');
fprintf('|----|-----------|-----------|-----------|\n');
fprintf('| %2d | %.6e | %.6f | %.6f |\n', iter, nor_grad, x_opt(1), x_opt(2));