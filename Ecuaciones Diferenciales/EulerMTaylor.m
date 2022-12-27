%*********Discretización de Euler *************
%*********       y'=f(t,y)        *************
clear all;
clc;

t0=0; %tiempo inicial
y0=1; %condición inicial y(t0)=y0
h=0.1; % paso
tmax=0.4; %tiempo de maximo

% comparación con la solución
t=t0:h:tmax;
y=1/2+1/2*exp(2*t);%solución exacta

ye(1)=y0;
ym(1)=y0;
yt(1)=y0;

for n=1:1:1+(tmax-t0)/h

    fe=2*y(n)-1;
    ye(n+1)=ye(n)+fe*h; %Euler Formula 
   
    fm=2*y(n)-1;
    fm1=2*(ym(n)+h*fm)-1;
    ym(n+1)=ym(n)+h*((fm+fm1)/2); %Euler Formula Mejorada

    ft=2*y(n)-1;
    ftd=2;
    yt(n+1)=yt(n)+h*(ft+ftd*(h/2)); %taylor con 3 terminos
end
plot(t,y,'b',t,ye,'r--',t,ym,'g:',t,yt);
grid on;

ee=norm(y-ye)
em=norm(y-ym)
et=norm(y-yt)

%Funciones ejemplo
% function f=pepito(t,y)
% f=1-t+4*y;
% end