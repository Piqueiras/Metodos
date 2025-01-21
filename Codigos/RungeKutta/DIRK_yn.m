function [x_vals, y_vals] = DIRK_yn(fid, f, a, eta, h, N, A, b, c, sol)
    % RUNGE KUTTA DIAGONALMENTE IMPLÍCITO GENERAL CON EL MÉTODO DE LOS yni

    % Entradas:
    %   f   - función que define el sistema de EDO
    %   a   - punto inicial
    %   eta - condición inicial
    %   h   - tamaño de paso
    %   N   - número de pasos 
    %       AVISO: No se implementa un b como punto final. 
    %       Hacer N = ceil((b-a)/N)
    %   A   - matriz de coeficientes de Butcher (tabla)
    %   b   - vector de pesos del método
    %   c   - vector de nodos del método
    %   sol - solución exacta (opcional)
    %
    % Salidas:
    %   x_vals - vector de los tiempos
    %   y_vals - solución numérica en cada tiempo

    % Inicialización de los vectores de tiempo y solución
    x_vals = a + (0:N) * h; 
    % Vector con todos los valores de x (a, a+h, a+2h...)
    y_vals = zeros(length(eta), N + 1);
    % Se toma la longitud de eta para considerar la dimensión
    % Se guardan entonces un vector de tamaño eta desde y0 hasta yN (N+1)
    y_vals(:,1) = eta;
    % Primer valor
    

    if nargin == 10
        err_vals = abs(y_vals(:,1)-sol(x_vals(1)));
    end

    % Imprimir cabecera
    if fid > 0
        if nargin == 10
            escribe_cabecera(fid, x_vals(1), y_vals(:,1), err_vals(:,1));
            escribe_paso(fid,0,x_vals(1),y_vals(:,1),err_vals(:,1));
        else
            escribe_cabecera(fid, x_vals(1), y_vals(:,1));
            escribe_paso(fid,0,x_vals(1),y_vals(:,1));
        end
    end

    s = length(c);
    yn = zeros(size(y_vals,1),s);
    % Iteración sobre los pasos
    for n = 1:N % Calcular en cada paso el yn
        for i = 1:s % Calcular en cada paso el yni
            sumaAK = 0;
            for j = 1:(i-1)
                sumaAK = sumaAK + A(i,j) * f(x_vals(n) + h*c(j),yn(:,j));
            end
            % sumaAK = sum(aij*f(xnj,ynj) hasta i-1
            % yni = f(xni, yn + h*sumaAK + h*(aii*f(xni,yni))
            % Se hace sistema con z siendo los yni
            yn(:,i) = fsolve(@(z) z - y_vals(:,n) - h*sumaAK - h*A(i,i)*f(x_vals(n)+h*c(i),z), y_vals(:,n), optimoptions(@fsolve,'Display','off','TolFun',1e-10,'TolX',1e-10));
        end

        % yn+1 = yn + h*Σ(bi*f(xni,yni))
        sumaBK = 0;
        for i = 1:s
            sumaBK = sumaBK + b(i) * f(x_vals(n) + h*c(i), yn(:,i));
        end
        y_vals(:,n+1) = y_vals(:,n) + h * sumaBK;

        if fid > 0
            if nargin == 10
                err_vals(:,n+1) = abs(y_vals(:,n+1) - sol(x_vals(n+1)));
                % Escribir el resultado con el error
                escribe_paso(fid, n, x_vals(n+1), y_vals(:,n+1), err_vals(:,n+1));
            else
                % Escribir el resultado sin error
                escribe_paso(fid, n, x_vals(n+1), y_vals(:,n+1));
            end
        end
        
    end
end