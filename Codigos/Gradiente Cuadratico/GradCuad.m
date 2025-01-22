function [x, nor_grad, index] = GradCuad(A, b, x0, itmax, tol)
    % GradCuad - Método de gradiente con paso óptimo
    %
    % Entradas:
    %   A      - Matriz simétrica definida positiva (n x n)
    %   b      - Vector (n x 1)
    %   x0     - Vector inicial (n x 1)
    %   itmax  - Máximo número de iteraciones
    %   tol    - Tolerancia para el test de parada
    %
    % Salidas:
    %   x         - Último iterante calculado
    %   nor_grad  - Norma del gradiente en el último iterante
    %   index     - Iteración en la que converge, -1 si no converge

    % Inicialización
    x = x0;
    gk = A * x - b; % Gradiente inicial
    nor_grad = norm(gk);
    index = -1; % Inicialmente asumimos que no converge
    
    % Iteraciones
    for k = 1:itmax
        if nor_grad <= tol
            index = k - 1; % Guardamos la iteración de convergencia
            break;
        end
        
        % Cálculo del paso óptimo
        alpha_k = (gk' * gk) / (gk' * A * gk);
        
        % Actualización del iterante
        x = x - alpha_k * gk;
        
        % Actualización del gradiente
        gk = A * x - b;
        
        % Actualización de la norma del gradiente
        nor_grad = norm(gk);
    end
end
