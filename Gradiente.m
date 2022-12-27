% Gradiente

% x^2-24*x+y^2-10*y
    % x = 8 , y = 7
    % x = 12 , y = 5

syms x y f  
f = input('Ingrese la funcion f(x,y) \n');
f = sym(f);
x0 = input('Ingrese el valor inicial para X \n');
y0 = input('Ingrese el valor inicial para Y\n');

dfdx = diff(f,x);
dfdy = diff(f,y);
d2fdx2 = diff(dfdx,x);
d2fdy2 = diff(dfdy,y);
d2fdxdy = diff(dfdx,y);
d2fdydx = diff(dfdy,x);

gradf = [dfdx; dfdy];   % Gradiente de la funcion
H = [d2fdx2 d2fdxdy; d2fdydx d2fdy2];   % Hessiano de la funcion
H1 = inv(H); %inversa de H

k = 0; %No. de intersecciones
e = 1; %Error
xi = x0;
yi = y0;

while e>= 1E-5 %0.00001
    x = xi;
    y = yi;
    gradf_ = eval(gradf);
    H1_ = eval(H1); %encuentra alpha, magnitud del paso
    M = [xi;yi]-H1_*gradf_; %para encontrar los valores del siguiente punto, o sea x1
    xi1 = M(1);
    yi1 = M(2);
    e = abs(sqrt((xi1-xi)^2 + (yi1-yi)^2)); %encuentra error
    xi = xi1; %nuevos puntos
    yi = yi1;
    k = k+1; %cuenta las iteraciones
    if k>100
        error('No se encontro solucion, k > 100 iteraciones')
    end
end

mensage = 'Se encontro una solucion';

x = xi;
y = yi;

detH = eval(det(H));
d2f = eval(d2fdx2);

%de acuerdo al determinante va a encontrar:

if detH < 0
    c = 'Punto silla';
else
    if d2f > 0
        c = 'Minimo';
    else
        c = 'Maximo';
    end
end
fprintf('Nuevos puntos de iteracion')
fprintf('X = %8.5f\nY = %8.5f\n',[xi,yi])
disp(' ')
fprintf('Numero de iteraciones necesarias:')
disp(k)
disp(' ')
fprintf('Se encontro un: ')
disp(c)
disp(' ')
