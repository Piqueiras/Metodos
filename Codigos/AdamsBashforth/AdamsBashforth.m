function [xn, yn] = AdamsBashforth(fid, f, a, eta, h, N, sol)
    format long;
    
    % Inicialización
    xn = a;
    yn = eta;
    Y = zeros(length(eta), N+1); % Almacena todas las soluciones
    X = a:h:a+h*N;               % Almacena los valores de x
    Y(:,1) = yn;                 % Condición inicial

    % Registro de resultados iniciales usando phi si existe
    if nargin == 7
        errn = Y(:,1) - sol(X(1));
        errmax = norm(errn);
        escribe_cabecera(fid, X(1), Y(:,1), errn);
        escribe_paso(fid, 0, X(1), Y(:,1), errn);
    else
        escribe_cabecera(fid, X(1), Y(:,1));
        escribe_paso(fid, 0, X(1), Y(:,1));
    end

    % Usar el método de Heun para los primeros dos valores
    % Esto establece y0, y1 y y2 necesarios para el método de tres pasos
    for i = 1:2
        k1 = f(X(i), Y(:,i));
        k2 = f(X(i) + h, Y(:,i) + h * k1);
        
        Y(:,i+1) = Y(:,i) + (h/2) * (k1 + k2);
        xn = X(i+1);

        if nargin == 7
            errn = Y(:,i+1) - sol(xn);
            errmax = max(errmax, norm(errn));
            escribe_paso(fid, i, xn, Y(:,i+1), errn);
        else
            escribe_paso(fid, i, xn, Y(:,i+1));
        end
    end

    fn = f(X(2),Y(:,2));
    fn_1 = f(X(1),Y(:,1));
    
    

    % Método de Adams-Bashforth de 3 pasos
    for n = 3:N
        % Evaluar f en los tres puntos previos
        fn_2 = fn_1;
        fn_1 = fn;
        fn = f(X(n), Y(:,n));
        
        % Adams-Bashforth de 3 pasos
        Y(:,n+1) = Y(:,n) + (h/12) * (23 * fn - 16 * fn_1 + 5 * fn_2);
        xn = X(n+1);  % Actualizar xn al valor correspondiente
        
        % Calcular error si se proporciona phi
        if nargin == 7
            errn = Y(:,n+1) - sol(xn);
            errmax = max(errmax, norm(errn));
            escribe_paso(fid, n, xn, Y(:,n+1), errn);
        else
            escribe_paso(fid, n, xn, Y(:,n+1));
        end
    end

    % Mostrar error máximo al final
    if nargin == 7
        fprintf(fid, "\nError máximo = %10.3e\n", errmax);
    end
    
    yn = Y(:,N+1); % Devolver la matriz completa de resultados
end
