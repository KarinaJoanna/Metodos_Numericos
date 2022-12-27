%****Discretización de Euler ******
%***       y'=f(t,y)        *****
clear all;
clc;

t0=0; %tiempo inicial
y0=1; %condición inicial y(t0)=y0
h=0.1; % paso
tmax=0.4; %tiempo de maximo

% comparación con la solución
t=t0:h:tmax;
%y=((19)/(16))*exp(4*t)+(1/4)*t-(3/(16)); %solución exacta

ye(1)=y0;
ym(1)=y0;
yt(1)=y0;
yr(1)=y0;
ymm(1)=y0;

for n=1:1:(tmax-t0)/h
    ye(n+1)=ye(n)+fn(t(n),ye(n))*h; %Euler Formula 
    
    ym1=ym(n)+h*fn(t(n),ym(n));
    ym(n+1)=ym(n)+h*(fn(t(n),ym(n))+fn(t(n+1),ym1))/2; %Euler Formula Mejorada

    yt(n+1)=yt(n)+h*(fn(t(n),yt(n))+fdn(t(n),yt(n))*(h/2)); %taylor con 3 terminos

    ymm(n+1)=ymm(n)+h*(fn(t(n)+(h/2),ymm(n)+((h/2)*fn(t(n),ymm(n))))); %Euler Modificado

    kn1=fn(t(n),yr(n));
    kn2=fn(t(n)+h/2,yr(n)+(kn1*h)/2);
    kn3=fn(t(n)+h/2,yr(n)+(kn2*h)/2);
    kn4=fn(t(n+1),yr(n)+(kn3*h));
    yr(n+1)=yr(n)+(h/6)*(kn1+2*kn2+2*kn3+kn4);%Runge-Kutta
end

plot(t,ye,'r--',t,ym,'k--',t,yt,'g:',t,yr,'m-.');
legend({'$ \ y_{E} \ $','$ \ y_{Em} \ $',' $ \ y_{T} \ $',' $ \ y_{R} \ $'},'interpreter','latex','FontSize',30,'Location','North','Orientation','Horizontal','EdgeColor',[0.93 0.93 0.93]);
% plot(t,y,'b',t,ye,'r--',t,ym,'k--',t,yt,'g:',t,yr,'m-.');
% legend({'$\ \phi (t) \ $','$ \ y_{E} \ $','$ \ y_{Em} \ $',' $ \ y_{T} \ $',' $ \ y_{R} \ $'},'interpreter','latex','FontSize',30,'Location','North','Orientation','Horizontal','EdgeColor',[0.93 0.93 0.93]);
grid on;

ee=norm(y-ye)
em=norm(y-ym)
et=norm(y-yt)
er=norm(y-yr)

% Funciones
function f=fn(t,y)
f=2*y-1;
end

function fd=fdn(t,y)
fd=2*t+2*y*fn(t,y);
end