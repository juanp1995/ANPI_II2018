#!/usr/bin/octave-cli

## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## II Semestre 2018
## Examen Final

## PROBLEMA 1

## NOMBRE: Juan Pablo Brenes Coto
## CARNE: 2014088353

1;

global m = 0.1;   ## Masa de la partícula
global b = 0.05;  ## Coeficiente de atenuación
global k = 1;     ## Constante de Hook

## ################## 
## ## Problema 1.1 ##
## ################## 
## Fuerza aplicada en la partícula
global F=@(x,v,t) -b*v -k*x;  ## <<< Ponga aquí su solución


## Resuelva el sistema atenuado masa resorte usando Euler
## tn el último instante de tiempo
## Dt paso temporal
function [t,x]=eulersys(tn,Dt)

  global m b k F;
  
  t=0:Dt:tn; ## Intervalo de simulación

  ## Pre-reserve la memoria utilizada.
  x=zeros(size(t));
  v=zeros(size(t));

  ## Condiciones iniciales
  x(1)=-1;

  ## ################## 
  ## ## Problema 1.2 ##
  ## ################## 

  ## Resuelva el sistema de ecuaciones con Euler

  ## >>> Ponga aquí su solución <<<
  for i = 1:length(t)-1
  	a = F(x(i), v(i), t(i)) / m; ## aceleración
  	
  	x(i+1) = x(i) + v(i)*Dt;
  	v(i+1) = v(i) + a*Dt;
  endfor
endfunction
  
figure(1,"name","Euler");
hold off;
[t,x]=eulersys(10,0.05);
plot(t,x,"r;\\Delta t=0.05;");

hold on;
[t,x]=eulersys(10,0.01);
plot(t,x,"m;\\Delta t=0.01;");

[t,x]=eulersys(10,0.001);
plot(t,x,"b;\\Delta t=0.001;");

xlabel("t");
ylabel("x(t)");
axis([0,10,-2,2]);
grid on;

## Resuelva el sistema de ecuaciones con Runge-Kutta 4to orden
## tn Último instante de tiempo
## Dt Paso temporal (delta t)
function [t,x] = rksys(tn,Dt)
  global m b k F;

  t=0:Dt:tn; ## Intervalo de simulación

  ## Pre-reserve la memoria utilizada.
  x=zeros(size(t));
  v=zeros(size(t));

  ## Condiciones iniciales
  x(1)=-1;

  ## ################## 
  ## ## Problema 1.4 ##
  ## ################## 

  ## >>> Ponga aquí su solución <<<
  h = @(x,v,t) -(b*v)/m -(k*x)/m;			## aceleración
  g = @(x,v,t) -F(x,v,t)/b - (k*x)/b;	## velocidad
  
  for i=1:length(t)-1
  	k1 = g(x(i), v(i), t(i))*Dt;
  	l1 = h(x(i), v(i), t(i))*Dt;
  	k2 = g(x(i)+k1/2, v(i)+l1/2, t(i)+Dt/2)*Dt;
  	l2 = h(x(i)+k1/2, v(i)+l1/2, t(i)+Dt/2)*Dt;
  	k3 = g(x(i)+k2/2, v(i)+l2/2, t(i)+Dt/2)*Dt;
  	l3 = h(x(i)+k2/2, v(i)+l2/2, t(i)+Dt/2)*Dt;
  	k4 = g(x(i)+k3, v(i)+l3, t(i)+Dt)*Dt;
  	l4 = h(x(i)+k3, v(i)+l3, t(i)+Dt)*Dt;
  	
  	x(i+1) = x(i) + (k1 + 2*k2 + 2*k3 + k4)*1/6;
  	v(i+1) = v(i) + (l1 + 2*l2 + 2*l3 + l4)*1/6;
  endfor

endfunction

figure(2,"name","RK");
hold off;
[t,x]=rksys(10,0.05);
plot(t,x,"r;\\Delta t=0.05;");

hold on;
[t,x]=rksys(10,0.01);
plot(t,x,"m;\\Delta t=0.01;");

[t,x]=rksys(10,0.001);
plot(t,x,"b;\\Delta t=0.001;");

xlabel("t");
ylabel("x(t)");
axis([0,10,-2,2]);
grid on;


