clear all;

tol = 1e-6;       % Tolerancia
itmax = 1000;     % Número máximo de iteraciones
fid = 1;

% Bucle para probar con cada valor de n
for n = [5, 10, 100]
    % Construcción de la matriz A
    A = 4 * eye(n) + diag(ones(n-1, 1), 1) + diag(ones(n-1, 1), -1);
    
    % Construcción del vector b
    b = zeros(n, 1);
    for i = 1:n-1
        b(i) = 2 * (-1)^(i+1);
    end
    b(1) = 3;             % Primer elemento
    b(end) = 3 * (-1)^(n+1);  % Último elemento

    % Iterante inicial
    x0 = zeros(n, 1);

    J = @(x) 0.5 * x' * A * x - b' * x;
    dJdx = @(x) A * x - b;

    % Llamada al método de gradiente
    [x, nor_grad, index] = GradWolfe(J, dJdx, x0, itmax, tol, 0.2, 0.4, fid);

    % Resultados
    fprintf('Resultados para n = %d:\n', n);
    fprintf('  Iteraciones hasta la convergencia: %d\n', index);
    fprintf('  Norma del gradiente en la convergencia: %.6e\n', nor_grad);
    escribe_v("Iterante final",x);
end
