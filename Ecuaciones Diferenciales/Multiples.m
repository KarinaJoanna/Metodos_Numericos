%*********Discretización de Euler *************
%*********       y'=f(t,y)        *************
clear all;
clc;

t0=0; %tiempo inicial
y0=1; %condición inicial y(t0)=y0
h=0.1; % paso
tmax=6; %tiempo de maximo

% comparación con la solución
t=t0:h:tmax;
y=((19)/(16))*exp(4*t)+(1/4)*t-(3/(16)); %solución exacta

ye(1)=y0;
yem(1)=y0;
yt(1)=y0;
yr(1)=y0;
ya(1)=y0;
ym(1)=y0;
yca(1)=y0;

for n=1:1:1+(tmax-t0)/h

    ye(n+1)=ye(n)+fn(t(n),ye(n))*h; %Euler Formula 
    
    yem1=yem(n)+h*fn(t(n),yem(n));
    yem(n+1)=yem(n)+h*(fn(t(n),yem(n))+fn(t(n+1),yem1))/2; %Euler Formula Mejorada


    yt(n+1)=yt(n)+h*(fn(t(n),yt(n))+fdn(t(n),yt(n))*(h/2)); %taylor con 3 terminos

    kn1=fn(t(n),yr(n));
    kn2=fn(t(n)+h/2,yr(n)+(kn1*h)/2);
    kn3=fn(t(n)+h/2,yr(n)+(kn2*h)/2);
    kn4=fn(t(n+1),yr(n)+(kn3*h));
    yr(n+1)=yr(n)+(h/6)*(kn1+2*kn2+2*kn3+kn4);%Runge-Kutta

    if n<=3
        ya(n+1)=yr(n+1);
        ym(n+1)=yr(n+1);
        yca(n+1)=ya(n+1);
    else
        ya(n+1)=ya(n)+(h/24)*(55*fn(t(n),ya(n))-59*fn(t(n-1),ya(n-1))+37*fn(t(n-2),ya(n-2))-9*fn(t(n-3),ya(n-3))); %Fórmula predictora de Adams-Bashforth
        ym(n+1)=ym(n-3)+((4*h)/3)*(2*fn(t(n),ym(n))-fn(t(n-1),ym(n-1))+2*fn(t(n-2),ym(n-2))); %Fórmula de Milne
        yca(n+1)=ya(n)+(h/24)*(9*fn(t(n+1),ya(n+1))+19*fn(t(n),ya(n))-5*fn(t(n-1),ya(n-1))+fn(t(n-2),ya(n-2))); %fórmula correctora de Adams-Moulton
    end

end
% plot(t,ye,'r--',t,yem,'k--',t,yt,'g:',t,yr,'m-.');
% legend({'$ \ y_{E} \ $','$ \ y_{Em} \ $',' $ \ y_{T} \ $',' $ \ y_{R} \ $'},'interpreter','latex','FontSize',30,'Location','North','Orientation','Horizontal','EdgeColor',[0.93 0.93 0.93]);


plot(t,y,'b',t,ye,'r--',t,yem,'k--',t,yt,'g:',t,yr,'m-.',t,ya,'-.c',t,ym,'--',t,yca,':');
legend({'$\ \phi (t) \ $','$ \ y_{E} \ $','$ \ y_{Em} \ $',' $ \ y_{T} \ $',' $ \ y_{R} \ $',' $ \ y_{A} \ $',' $ \ y_{M} \ $',' $ \ y_{cA} \ $'},'interpreter','latex','FontSize',24,'Location','North','Orientation','Horizontal','EdgeColor',[0.93 0.93 0.93]);
grid on;

ee=norm(y-ye)
eem=norm(y-yem)
et=norm(y-yt)
er=norm(y-yr)
ea=norm(y-ya)
em=norm(y-ym)
eca=norm(y-yca)

% Funciones
function f=fn(t,y)
f=-1000*y(n)+5000;
end

function fd=fdn(t,y)
fd=-1000;
end