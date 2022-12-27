%*Discretización de Euler ***
%*       y'=f(t,y)        ***
clear all;

t0=0; %tiempo inicial
y0=1; %condición inicial y(t0)=y0
h=0.1; % paso
tmax=6; %tiempo de maximo

y(1)=y0;
t(1)=t0;

for n=1:1:tmax/h
    f=-1000*y(n)+5000;
    y(n+1)=y(n)+f*h; %Euler Formula
    t(n+1)=n*h+t0;
end

plot(t,y);