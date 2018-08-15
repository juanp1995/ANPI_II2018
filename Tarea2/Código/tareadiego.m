#Instituto Tecnologico de Costa Rica
#Area Academica de Ingenieria en Computadores
#Analisis Numerico para Ingenieria
#Estudiante: Diego Alonso Granados Ly
#Tarea_2

%Graficos con OpenGL
graphics_toolkit('gnuplot');
%Aumentamos la precision de salida
output_precision(20);
%Entrada del sistema.
xi = input("Ingrese un valor para x: ")
%Declaracion de las variables
h0     = 1;
lambda = double(0.9);
xil  = 0;%xi-1
xip  = 0;%xi+1
fxil = 0;%f(xi-1)
fxi  = 0;%f(xi)
fxip = 0;%f(xi+1)
vector_h   = zeros(1,500);%Vector con los tamanos de paso
vector_da  = zeros(1,500);%Vector con la diferencia hacia atras
vector_dd  = zeros(1,500);%Vector con la diferencia hacia adelante
vector_c   = zeros(1,500);%Vector con la diferencia centrada
error_da   = zeros(1,500);%Vector con el error de la diferencia hacia adelante
error_dd   = zeros(1,500);%Vector con el error de la diferencia hacia atras
error_c    = zeros(1,500);%Vector con la diferencia centrada

#Ciclo Que define los valores de xi , xi+1 y xi-1
#ademas de sus evaluaciones en la funcion.
for i=1:500
  if (i==1)
    vector_h(i) = h0;               %Excepcion para i = 1 con  h0 = 1
  else
    vector_h(i) = lambda*vector_h(i-1); %Con i!=1 se aplica h(i)=lambda*h(i-1) 
  endif
  xil  = xi - vector_h(i);      %xi-1 = xi-h
  xip  = xi + vector_h(i);      %xi+1 = xi+h
  fxil = double(sin(e^(xil^2.0)));%f(xi-1)
  fxi  = double(sin(e^(xi^2.0)));  %f(xi)
  fxip = double(sin(e^(xip^2.0)));%f(xi+1)

  vector_da(i) = double((fxi-fxil)/(xi-xil));
  vector_dd(i) = double((fxip-fxi)/(xip-xi));
  vector_c (i) = double((fxip-fxil)/(2.0*vector_h(i)));

  if (i==1)
    error_da(i)  = 0;
    error_dd(i)  = 0;
    error_c (i)  = 0;
  else
    error_da(i)  = abs(((vector_da(i-1)-vector_da(i))/(vector_da(i)))*100);
    error_dd(i)  = abs(((vector_dd(i-1)-vector_dd(i))/(vector_dd(i)))*100);
    error_c (i)  = abs(((vector_c (i-1)-vector_c (i))/(vector_c(i)))*100);
  endif
  
endfor

#Definimos las propiedades de la grafica.
loglog(vector_h,error_c,'r',vector_h,error_da,'g',vector_h,error_dd,'b');
axis([10^-15,10^0,10^-13,10^01]);
xlabel("Diferencia h");
ylabel("Error");
title("Aproximacion por diferencias");
legend("Centrada","Delante","Detras");
grid on;