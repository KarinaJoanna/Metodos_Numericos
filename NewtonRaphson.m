%Damos de alta la variable simbólica X
syms x
%Introducimos la función,el punto de inicio,así como
%porcentaje de error
f=input('Introduzca la función f(x):');
pi=input('Introduzca el punto de inicio:');
err=input('Porcentaje de error:');
%Graficamos la función
ezplot(f)
grid on
%Calculamos la derivada de la función
d=diff(f);
d=inline(d);
f=inline(f);
ea=100;
j=0;
while ea>err
%Aproximamos la raiz con la fórmula correpondiente
xi=pi-(f(pi)/d(pi));
%Calculamos el porcentaje de error
ea=abs(((xi-pi)/xi)*100);
pi=xi;
j=j+1;
end
%Mostramos los resultados en pantalla (con 3 decimales)
fprintf('\nRaiz= %10.3f en %d Iteraciones',pi,j)