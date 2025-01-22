function [x, nor_grad, index] = GradWolfe(J, dJdx, x0, itmax, tol, rho, sigma, fid)
    % Entradas:
    %   J       - Funcional coste
    %   dJdx    - Gradiente del funcional
    %   x0      - Vector inicial (n x 1)
    %   itmax   - Máximo número de iteraciones
    %   tol     - Tolerancia para el test de parada
    %   rho & sigma
    %
    % Salidas:
    %   x         - Último iterante calculado
    %   nor_grad  - Norma del gradiente en el último iterante
    %   index     - Iteración en la que converge, -1 si no converge
    
    % Inicialización
    alphainit = 1;
    gamma = 2;
    x = x0;
    gk = dJdx(x); % Gradiente inicial
    nor_grad = norm(gk);
    dk = -gk;
    index = -1; % Inicialmente asumimos que no converge

    fprintf(fid, '\nIteración\t||grad||\t\tJ(x)\t\tAlpha\n');
    fprintf(fid, '--------------------------------------------------\n');
    
    % Iteraciones
    for k = 1:itmax
        if nor_grad <= tol
            index = k - 1; % Guardamos la iteración de convergencia
            break;
        end
        
        [alpha_k, ~, index] = WolfeBiseccion(J, dJdx, x, dk, rho, sigma, alphainit, gamma, fid);
        
        % Imprimir información de la iteración actual
        fprintf(fid, '%4d\t\t%.6e\t%.6e\t%.6e\n', k, nor_grad, J(x), alpha_k);
        
        % Actualización del punto
        x = x + alpha_k * dk;
        
        % Actualización del gradiente
        gk = dJdx(x);
        
        % Actualización de la dirección de descenso (Gradiente descendente estándar)
        dk = -gk;
        
        % Actualización de la norma del gradiente
        nor_grad = norm(gk);
    end
    
    if nor_grad > tol
        fprintf(fid, 'No se alcanzó la tolerancia en %d iteraciones.\n', itmax);
    else
        fprintf(fid, 'Convergencia alcanzada en %d iteraciones.\n', index);
    end
end
