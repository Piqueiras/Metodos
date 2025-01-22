function [alpha, i0, index] = WolfeBiseccion(f, g, x, d, ...
                                    rho, sigma, ...
                                    alphainit, gamma, ...
                                    fid)
                                
% Wolfe   Seleccion de paso para un metodo de descenso mediante la regla
%         inexacta llamada regla de Wolfe. La generacion de un nuevo punto
%         alpha entre alpha_p y alpha_g se realiza por biseccion.
%
% Parametros de entrada:
%     f         : Funcion coste (psoblemente multi-dimensional).
%     g         : Gradiente de la funcion coste.
%     x         : Punto en el que se encuentra el algoritmo de descenso.
%     d         : Direccion de descenso.
%     rho, sigma: Parametros que definen los criterios de la regla de 
%                 Wolfe. Se supone que 0 < rho < sigma < 1.
%                 Es aconsejable que   0 < rho < 1/2.
%     alphainit : Valor de alpha que se considera inicialmente.
%                 Deberia satisfacer el criterio de paso grande.
%     gamma     : Parametro empleado para seleccionar un paso inicial que 
%                 verifique el criterio de paso grande en caso de que 
%                 alphag no lo cumpla.
%                 Se trata de un factor de dilatacion. En consecuencia debe
%                 tomarse tal que gamma > 1.
%     fid       : Identificador de un fichero para la salida.
%
% Parametros de salida:
%     alpha     : Paso seleccionado por la regla de Wolfe para el metodo de
%                 descenso.
%     i0        : Numero de evaluaciones del funcinal coste.
%     index     : Indicador de convergencia: 1 => exito, <= => fracaso.
%                 Si index = 0  => No se ha podido inicializar.
%                 Si index = -1 => No se ha podido encontrar un alpha
%                 admisible.

% Inicializacion del algoritmo. 
% Seleccionamos un alpha que cumpla el criterio de paso grande y otro que
% cumpla el criterio de paso pequeno.

index  = 0; 
i0     = 0;   % Numero de evaluaciones del funcional coste.
itmax1 = 10;  % Numero maximo de iteraciones en la inicializacion. 
itmax2 = 50;  % Numero maximo de iteraciones en la busqueda. 

init   = false;
alphap = 0;
alpha  = alphainit;

j0  = f(x);
jp0 = g(x).'*d;

while ~init && ( i0 <= itmax1 )
    xalpha  = x + alpha*d;
    jalpha  = f(xalpha);
    jpalpha = g(xalpha).'*d;
    
    CPG = (jalpha > j0 + rho * jp0 * alpha); % Criterio de paso grande.
    CPP = (jpalpha < sigma * jp0) && (~CPG); % Criterio de paso pequeno.
    
    i0 = i0 + 1;
    
    if CPG                         % Alpha cumple CPG
        alphag   = alpha;
        init     = true;
    elseif CPP                     % Alpha cumple CPP
        alphap   = alpha;
        alpha    = alpha * gamma;
    else                           % Alpha cumple CPA
        index = 1;

        return;                    % Salimos con exito de la funcion
    end
end

if (i0 > itmax1) && (init == false)    % No se ha podido inicializar Wolfe.
    return;
end

% Algoritmo inicializado. Disponemos de un alphap que cumple el criterio
% de paso pequeno y un alphag que cumple el criterio de paso grande.
% Buscamos un paso admisible entre ambos.

while ( i0 <= itmax2 )
    
% Seleccionamos el nuevo alpha con el metodo de biseccion.

    alpha = (alphag + alphap)/2;
            
    xalpha  = x + alpha*d;
    jalpha  = f(xalpha);
    jpalpha = g(xalpha).'*d;

    CPG = (jalpha > j0 + rho * jp0 * alpha); % Criterio de paso grande.
    CPP = (jpalpha < sigma * jp0) && (~CPG); % Criterio de paso pequeno.

    i0 = i0 + 1;
    
    if CPG                         % Alpha cumple CPG
        alphag   = alpha; 
    elseif CPP                     % Alpha cumple CPP
        alphap   = alpha;
    else                           % Alpha cumple CPA
        index = 1;
        return;
    end
end

index = -1;
return;

end


