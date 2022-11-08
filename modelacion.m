% Equipo 3
% 
% ATENCION: ESTE SCRIPT FUNCIONA. 
% 
% NO MOVER NADA A ESTE DOCUMENTO A MENOS QUE SE MUESTRE QUE EL DOCUMENTO DE
% PRUEBAS FUNCIONE

%{

PARAMETROS: CAMBIAN SEGUN LOS DATOS INGRESADOS

v0 = 20; % Velocidad Inicial (m/s)
theta = 40; % Angulo entre la Velocidad inicial y la horizontal, en grados
dt = 0.05; % Tiempo entre puntos medidos 
coef = 0.45; % Coeficiente de Arrastre
ro = 1.05; % Densidad del aire (kg/m^3)
m = 23.168; % Masa (kg)
radio = 0.12; % Metros (m)

%}

function [x,y,xmax,ymax,t_total,e,constante] = modelacion(v0,theta,dt,coef,ro,m,radio)
    clc

    % CONSTANTES: NO SE MUEVEN NUNCA
    x0 = 0; % Posicion del objeto en el eje de X inicial (m)
    y0 = 20; % Posicion del objeto en el eje de Y inicial (m)
    g = 9.81; % Gravedad (m/s^2)
    t0 = 0; % Tiempo inicial (s)

    % OPERACIONES: Crear listas para modelar el movimiento
    a = pi()*radio^2; % Area (m^2)
    constante = coef*ro*a/(2*m); % Constante de Arrastre
    
    % Genera la posicion en el eje de X y Y (con resistencia)
    [t,x,y] = lista(v0,theta,t0,dt,x0,y0,constante,g);
    
    x;
    y;
    xmax = max(x);
    ymax = max(y);

    t_total = max(t);

    e = error(ymax,v0,theta,t0,x0,y0,constante,g);

    render(x,y)
end

function [t,x,y] = lista(v0,theta,t0,dt,x0,y0,constante,g)
    % Componentes de X y Y del vector de velocidad inicial
    v0x = v0*cos(theta*pi/180); % Componente de X de la Velocidad inicial (m/s)
    v0y = v0*sin(theta*pi/180); % Componente de Y de la Velocidad inicial (m/s)
    
    % El primer elemento de la lista siempre sera el tiempo inicial.
    t(1) = t0;
    
    % El primer elemento de la lista X y Y siempre seran sus 
    % respectivas posiciones iniciales.    
    x(1) = x0;
    y(1) = y0;

    % El primer elemento de la lista vx y vy siempre seran sus 
    % respectivas velocidades iniciales.  
    vx(1) = v0x;
    vy(1) = v0y;

    % El primer elemento de la lista ax y ay siempre seran sus 
    % respectivas aceleraciones iniciales tomando en cuenta la resistencia.
    ax(1) = -constante*v0*v0x;
    ay(1) = -g-constante*v0*v0y;
    
    % Índice
    i=1;
    
    % Método de Euler, se detiene cuando la coordenada en Y es menor a 0. 
    while y(i) > 0
        t(i+1) = t(i) + dt;
        
        x(i+1) = x(i) + vx(i)*dt;
        y(i+1) = y(i) + vy(i)*dt;
    
        vx(i+1) = vx(i) + ax(i)*dt;
        vy(i+1) = vy(i) + ay(i)*dt;
    
        ax(i+1) = -constante*sqrt(vx(i+1)^2+vy(i+1)^2)*vx(i+1);
        ay(i+1) = -g-constante*sqrt(vx(i+1)^2+vy(i+1)^2)*vy(i+1);
    
        i = i + 1;
    end 
        
    % Interpolación de (x,0), cuando el objeto hace contacto con la
    % horizontal
    n = length(y);
    y(n) = 0;
    x(n) = -y(n-1)/(vy(n-1)/vx(n-1)) + x(n-1);
    t(n) = t(n-1) + (x(n)-x(n-1))/vx(n-1);

end

function e = error(ymax,v0,theta,t0,x0,y0,constante,g)
    [T,X,Y] = lista(v0,theta,t0,0.00001,x0,y0,constante,g);
    h_real = max(Y); % dt = 0.00001
    e = 100*abs(h_real-ymax)/h_real;
end

function cono()
    % Cono
    R = 1000; % Radio base
    r = 500; % Radio altura
    
    H = 40; 
    
    ncs = 10;
    nNodes = 100;
    
    z=linspace(-H/2,H/2,ncs);
    m = -(R-r)/H;
    b = (R+r)/2;
    r_local = m*z+b;
    th = linspace(0,2*pi,nNodes);
    
    [th, r_local] = meshgrid(th, r_local);
    [X,Y] = pol2cart(th,r_local);
    
    Z = repmat(z',1,nNodes);
    
    surf(X,Y,Z)
end

function render(x,y)
    % Render en gráfica de los puntos considerando resistencia
    figure(1)
    plot(x,y,"o")
    grid
    
    % Leyenda de que significa el color de los puntos
    legend('Trayectoria con Resistencia','Location','southwest')
    
    % Formato de gráfica
    title('Trayectoria de Objeto Expulsado')
    xlabel('Distancia Horizontal (metros)') 
    ylabel('Distancia Vertical (metros)')
    ylim([0,max(y)+1])
    xlim([0,max(x)])
    
    figure(2)
    phi = 0:0.01:2*pi;
    
    for i = 10:5:25
        [t,x,y] = lista(i,40,0,0.05,0,20,0.00046131,9.81);
        max(x)
        p = polarplot(phi,max(x)*(sin(phi).^2+cos(phi).^2));
        color = [1-(i-10)/25 (i-10)/25 0];
        p.Color = color;
        p.LineWidth = 3;
        hold on
    end

    title("Mapa de Riesgo")
    %rticks([21.4,37.2,56.7,80.1]);
    %rticklabels({'V = 10','V = 15','V = 20','V = 25'});

    
    figure(3)
    temp = [1:1:length(x)];
    plot3(x,temp,y)
    hold on
    cono()
    title('Modelo en 3 Dimensiones de Movimiento')
    xlabel('Longitud (metros)') 
    ylabel('Latitud (metros)')
    zlabel('Altitud (metros)')
    axis equal
    grid

end