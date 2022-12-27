%*********Discretización de Euler *************
%*********       y'=f(t,y)        *************
clear all;
clc;

t0=0; %tiempo inicial
y0=[1;0]; %condición inicial y(t0)=y0
h=0.1; % paso
tmax=t0+2; %tiempo de maximo


% comparación con la solución
t=t0:h:tmax;

ye(:,1)=y0;
yem(:,1)=y0;
yt(:,1)=y0;
yr(:,1)=y0;
ya(:,1)=y0;
ym(:,1)=y0;
yca(:,1)=ya(:,1);


for n=1:1:(tmax-t0)/h

    ye(:,n+1)=ye(:,n)+fn(t(n),ye(:,n))*h; %Euler Formula 
    
    yem1=yem(:,n)+h*fn(t(n),yem(:,n));
    yem(:,n+1)=yem(:,n)+h*(fn(t(n),yem(:,n))+fn(t(n+1),yem1))/2; %Euler Formula Mejorada


    yt(:,n+1)=yt(:,n)+h*(fn(t(n),yt(:,n))+fdn(t(n),yt(:,n))*(h/2)); %taylor con 3 terminos

    kn1=fn(t(n),yr(:,n));
    kn2=fn(t(n)+h/2,yr(:,n)+(kn1*h)/2);
    kn3=fn(t(n)+h/2,yr(:,n)+(kn2*h)/2);
    kn4=fn(t(n+1),yr(:,n)+(kn3*h));
    yr(:,n+1)=yr(:,n)+(h/6)*(kn1+2*kn2+2*kn3+kn4);%Runge-Kutta

    if n<=3
        ya(:,n+1)=yr(:,n+1);
        ym(:,n+1)=yr(:,n+1);
        yca(:,n+1)=ya(:,n+1);
    else
        ya(:,n+1)=ya(:,n)+(h/24)*(55*fn(t(n),ya(:,n))-59*fn(t(n-1),ya(:,n-1))+37*fn(t(n-2),ya(:,n-2))-9*fn(t(n-3),ya(:,n-3))); %Fórmula predictora de Adams-Bashforth
        ym(:,n+1)=ym(:,n-3)+((4*h)/3)*(2*fn(t(n),ym(:,n))-fn(t(n-1),ym(:,n-1))+2*fn(t(n-2),ym(:,n-2))); %Fórmula de Milne
        yca(:,n+1)=ya(:,n)+(h/24)*(9*fn(t(n+1),ya(:,n+1))+19*fn(t(n),ya(:,n))-5*fn(t(n-1),ya(:,n-1))+fn(t(n-2),ya(:,n-2))); %fórmula correctora de Adams-Moulton
    end

end

plot(t,ye,'--',t,yem,'-',t,yt,':',t,yr,'-.',t,ya,'--',t,ym,'-',t,yca,':');
grid on;

% Funciones
function f=fn(t,y)
f1=y(2);
f2=y(3);
f3=-y(3)-y(2)+y(1)+exp(-t)-t^2+2*t+2;
f=[f1; f2; f3];
end

function fd=fdn(t,y)
fd1=y(3);
fd2=-y(3)-y(2)+y(1)+exp(-t)*cos(t)-t^2+2*t+2;
fd3=-y(1)-exp(-t)*(2*cos(t)+sin(t))-4*t*t^2;
fd=[fd1; fd2; fd3]
end