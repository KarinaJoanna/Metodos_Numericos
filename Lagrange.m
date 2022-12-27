% Lagrange
   % f = 9-x^2-y^2
   % g = x+y-3

   % f = x*y+2*x
   % g = 4*x+2*y - 60

i = 2;
syms x y f g l
f(x,y) = input('Ingrese la funcion f(x,y) \n'); % Funcion objetivo
f = sym(f);

g(x,y) = input('Ingrese la restriccion g(x,y)\n ');  % Funcion de restriccion
g = sym(g);

%encontrar las expresiones

fxdx = diff(f(x,y),x);
fydy = diff(f(x,y),y);

fx(x,y) = (-1 * fxdx )*(l * diff(g,x));
fy(x,y) = (-1 * fydy )*(l * diff(g,y));

solve(fx(x,y),l);
solve(fy(x,y),l);
solve(g(x,y),x);

L(x,y,l) = f + l*g(x,y);%lagrange
%L(x,y,l) = f(x,y) + l*g(x,y);

dgdx = diff(g(x,y),x);
dgdy = diff(g(x,y),y);
dLdx = diff(L(x,y,l),x);
d2Ldx2 = diff(dLdx,x);
dLdy = diff(L(x,y,l),y);
d2Ldy2 = diff(dLdy,y);
d2Ldxdy = diff(dLdx,y);
d2Ldydx = diff(dLdy,x);

%matriz hessiano orlando
HO(x,y,l) = [0 dgdx dgdy;
      dgdx d2Ldx2 d2Ldydx;
      dgdy d2Ldxdy d2Ldy2];

HOD = det(HO(x,y,l));
max = HOD * (-1)^(i+1);

if HOD < 0 
    c = 'Minimo';
else
    if max < 0 
        c = 'Maximo';
    else
        c = 'Punto silla';
    end
end

%determina si es un maximo o un minimo
fprintf('Se encontro un: ')
disp(c)
disp(' ')