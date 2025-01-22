function [x, nor_grad, index] = GradConjCuad(A, b, x0, itmax, tol)
    % GradConjCuad - Método de gradiente conjugado con paso óptimo
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
    dk = -gk;       % Primera dirección de descenso
    nor_grad = norm(gk);
    index = -1;     % Inicialmente asumimos que no converge

    % Iteraciones
    for k = 1:itmax
        Adk = A * dk;
        alpha = (gk' * gk) / (dk' * Adk); % Cálculo del paso óptimo
        x = x + alpha * dk; % Actualización de x
        gk_new = gk + alpha * Adk; % Se necesita el nuevo y anterior

        nor_grad = norm(gk_new);
        if nor_grad < tol
            index = k; % Convergencia
            break;
        end

        beta = (gk_new' * gk_new) / (gk' * gk); % Cálculo de beta
        dk = -gk_new + beta * dk; % Nueva dirección de descenso
        gk = gk_new;
    end

    if nor_grad >= tol
        index = -1; % No convergió en itmax iteraciones
    end
end
