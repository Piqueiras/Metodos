function [xN, yN] = RK4(fid, f, a, eta, h, N, sol)
    % fid - identificador de fichero para la salida (1 para pantalla)
    % f   - función que define el sistema de EDO
    % a   - punto inicial
    % eta - condición inicial
    % h   - tamaño de paso
    % N   - número de pasos
    % sol - solución exacta (opcional)
    
    % Inicialización de los vectores
    x_vals = a + (0:N) * h;
    y_vals = zeros(length(eta), N + 1);  % Tamaño basado en eta (puede ser vector)
    y_vals(:,1) = eta;

    % Inicialización de errores si se proporciona la solución exacta
    if nargin == 7
        err_vals = zeros(size(y_vals));
        for k = 1:length(sol)
            err_vals(k,1) = abs(y_vals(k,1) - sol{k}(x_vals(1)));
        end
    end
    
    % Imprimir cabecera
    if fid > 0
        if nargin == 7
            escribe_cabecera(fid, x_vals(1), y_vals(:,1), err_vals(:,1));
            escribe_paso(fid,0,x_vals(1),y_vals(:,1),err_vals(:,1));
        else
            escribe_cabecera(fid, x_vals(1), y_vals(:,1));
            escribe_paso(fid,0,x_vals(1),y_vals(:,1));
        end
    end
    
    
    % Iteración sobre los pasos
    for n = 1:N
        % Cálculo de las pendientes k1, k2, k3, k4
        k1 = f(x_vals(n), y_vals(:,n));
        k2 = f(x_vals(n) + h/2, y_vals(:,n) + (h/2) * k1);
        k3 = f(x_vals(n) + h/2, y_vals(:,n) + (h/2) * k2);
        k4 = f(x_vals(n) + h, y_vals(:,n) + h * k3);

        % Actualización del valor de y en el siguiente paso
        y_vals(:,n+1) = y_vals(:,n) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);

        % Si se proporciona la solución exacta, calculamos el error
        if fid > 0
            if nargin == 7
                for k = 1:length(sol)
                    err_vals(k,n+1) = abs(y_vals(k,n+1) - sol{k}(x_vals(n+1)));
                end
                % Escribir el resultado con el error
                escribe_paso(fid, n, x_vals(n+1), y_vals(:,n+1), err_vals(:,n+1));
            else
                % Escribir el resultado sin error
                escribe_paso(fid, n, x_vals(n+1), y_vals(:,n+1));
            end
        end
    end
  

    % Valores finales
    xN = x_vals(end);
    yN = y_vals(:,end);
end
