## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## I Semestre 2018
## Examen Final

## PROBLEMA 1

## NOMBRE: Jairo Ortega Calderón
## CARNE: 2014043224

2;

## Cargar los datos
X=load("-ascii","regresion.dat");
output_precision(8) ## Cambio la precisión del número
global h; ## Variable global para el 

####################################################
## Problema 1.1                                   ##
## Grafique los puntos bidimensionales            ##
####################################################
figure(1);
plot(X(1,1:end),X(2,1:end),'.' );
xlabel("X")
ylabel("Y")
grid on;
hold off;

####################################################
## Problema 1.2                                   ##
## Implemente la función de error                 ##
####################################################
function val=f(abc,X)
  ## abc: vector columna [a,b,c]' con los parámetros de la función cuadrática
  ## X:   datos para evaluar la función, un dato por columna
  
  ## Ponga su código aquí:
  iSum=1;
  tempVal=0;
  mSum=columns(X);
  while(iSum<=mSum)
    tempVal+=(X(2,iSum) - (abc(1)*(X(1,iSum)^2) + abc(2)*X(1,iSum) + abc(3)))^2;
    iSum+=1;
  endwhile
  val=tempVal;

endfunction


####################################################
## Problema 1.3                                   ##
## Implemente el gradiente de la función de error ##
####################################################
function val=gf(abc,X)
  ## abc: vector columna [a,b,c]' con los parámetros de la función cuadrática
  ## X:   datos para evaluar la función, un dato por columna
  
  ## Use diferenciación NUMERICA para calcular el gradiente de f:
  iSum=1;
  tempVal=[0,0,0];
  mSum=columns(X);
  global h;
  while(iSum<=mSum)
    tempVal+=diferenciaCentrada(abc,X(1,iSum), X(2,iSum), h);
    iSum+=1;
    
  endwhile
  val=tempVal';

endfunction

####################################################
## Problema 1.4                                   ##
## Descenso de gradiente                          ##
####################################################
function [ABC,err]=optimice(f,gf,X,lambda,tol,abc0=[0,0,0]')
  ## f      es el handler de la función a optimizar
  ## gf     es el handler que calcula el gradiente de f
  ## X      es la matriz de datos 
  ## lambda es el tamaño de paso del descenso de gradiente
  ## tol    es el umbral de tolerancia para determinar convergencia
  ## abc0   es un vector [a0,b0,c0] especificando el punto inicial de
  ##        la optimización
  ## ABC    es una matrix de n x 3, donde cada fila corresponde a un
  ##        paso en el proceso de optimización.  Es decir, ABC(:,1)
  ##        corresponde siempre a abc0, y ABC(:,rows(ABC)) corresponde
  ##        a los parámetros óptimos.
  ## err    es el vector conteniendo los errores en cada paso

  if ( (rows(abc0)!=3) || columns(abc0)!=1 )
    error("Vector inicial no tiene forma 3x1");
  endif;
 
  ## Ponga su código aquí:
  errorAprox=[];
  ABCAprox=[];
  ea=errorPorcentualAproximado(f(abc0', X),0);
  aproxAct=abc0';
  global h;
  h=lambda;
  n=1;
  while(abs(ea)>tol)
  ABCAprox=[ABCAprox;aproxAct];
  errorAprox=[errorAprox;abs(ea)];
  aproxAct=ABCAprox(rows(ABCAprox),1:end)-lambda*gf(ABCAprox(rows(ABCAprox),1:end),X)';
  ea=errorPorcentualAproximado(f(aproxAct,X),f(ABCAprox(rows(ABCAprox),1:end),X));
  h=h*0.5;
  n=n+1;
  endwhile
  err=errorAprox;
  ABC=ABCAprox;
endfunction

## Llame al optimizador con la interfaz anterior

lambda=0.001;  # Ajuste esto
tol=0.00000001; # Ajuste esto
[ABC,err]=optimice(@f,@gf,X,lambda,tol,[0,1,0]');

####################################################
## Problema 1.5                                   ##
## Imprima el conjunto óptimo de parámetros       ##
####################################################
ABCOptimo=ABC(rows(ABC),1:end);
a=ABCOptimo(1)
b=ABCOptimo(2)
c=ABCOptimo(3)
####################################################
## Problema 1.6                                   ##
## Muestre el error en función de las iteraciones ##
####################################################

figure(2);
plot(1:rows(err),err(1:end,1),'.' );
title('Error Porcentual Aproximado vs Número de iteración')
xlabel('Número de iteración')
ylabel('Error Porcentual Aproximado')
grid on;
hold off;

####################################################
## Problema 1.7                                   ##
## Muestre las curvas inicial, intermedias y      ##
## final ajustadas a los datos                    ##
####################################################

cont=2;
figure(3);
axis([-2,2, -4,4]);
plot(X(1,1:end),X(2,1:end),'xb' );
hold on;
Xabc=(-2:0.1:2);
Yabc0=(ABC(1,1:end)(1)*(Xabc.^2))+(ABC(1,1:end)(2)*Xabc)+(ABC(1,1:end)(3));
graficoBase=plot(Xabc,Yabc0,'k');
set(graficoBase, 'LineWidth',2);
hold on;
while(cont<(rows(ABC)-2))
hold on;
  YabcInt=(ABC(cont,1:end)(1)*(Xabc.^2))+(ABC(cont,1:end)(2)*Xabc)+(ABC(cont,1:end)(3));
  plot(Xabc,YabcInt,'c');
  cont+=1;
endwhile
hold on;
YabcOpt=(ABC(rows(ABC),1:end)(1)*(Xabc.^2))+(ABC(rows(ABC),1:end)(2)*Xabc)+(ABC(rows(ABC),1:end)(3));
graficoOpt=plot(Xabc,YabcOpt,'r');
set(graficoOpt, 'LineWidth',2);
hold on;
title('Resultado Final')
grid on;
xlabel('x')
ylabel('y')
hold off;
##################################################################################