#!/usr/bin/octave-cli

## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## II Semestre 2018
## Examen Final

## PROBLEMA 2

## NOMBRE: Juan Pablo Brenes Coto
## CARNE: 2014088353
2;


## Construya algunos datos 2D para el problema
## N: Número de datos.  Cada fila de la matriz tendrá un punto.
##    La primera columna tendrá la coordenada x y la segunda columna
##    la coordenada y.
function points = createData(N)
  astep = 360/N;
  
  angles = (0:astep:360-astep)';
  anoise = rand(size(angles))*astep/4;
  rnoise = rand(size(angles))*1.5;
  
  radii  = 2+rnoise;
  angles = deg2rad(angles + anoise);
  
  points = [radii.*cos(angles) radii.*sin(angles)];
endfunction

## Calcule las segundas derivadas
## x: posiciones x de las muestras
## f: valores de la función en cada x
## retorne fpp con los valores de la segunda derivada en cada posición t
function fpp=findDerivs(t,f)
  assert(size(f)==size(t));

  N=length(t)-1; # Número de subintervalos
  
  ## Arme el sistema de ecuaciones
  M=eye(N+1,N+1);
  fpp=zeros(N+1,1);
  b  =zeros(N+1,1);

  ## ################## 
  ## ## Problema 2.4 ##
  ## ################## 

  ## >>> Ponga su solución aquí <<<
    for i=1:N+1
		j=i;
		if(i==1) ## Punto i-1 debe ser N+1
			M(i,j) = 2*(t(i+1)-t(N+1));
			M(i,j+1) = t(i+1)-t(i);
			M(i,N+1) = t(i)-t(N+1);
			b(i) = 6*( (f(i+1)-f(i))/(t(i+1)-t(i)) ) - 6*( (f(i)-f(N+1))/(t(i)-t(N+1)) );
			
		elseif(i==(N+1)) ## Punto i+1 debe ser 1 (i-N)
			M(i,i-N) = t(i-N)-t(i);
			M(i,j) = 2*(t(i-N)-t(N));
			M(i,N) = t(i)-t(i-1);
			b(i) = 6*( (f(i-N)-f(i))/(t(i-N)-t(i)) ) - 6*( (f(i)-f(i-1))/(t(i)-t(i-1)) );
			
		else
			M(i,j-1) = t(i)-t(i-1);
			M(i,j) = 2*(t(i+1)-t(i-1));
			M(i,j+1) = t(i+1)-t(i);
			b(i) = 6*( (f(i+1)-f(i))/(t(i+1)-t(i)) ) - 6*( (f(i)-f(i-1))/(t(i)-t(i-1)) );
		endif		
  endfor

  ## Resuelva el sistema
  fpp = M\b; 

endfunction

## Interpole los valores fi(xi) usando los puntos x y sus valores f(x)
## t: valores de soporte conocidos
## f: valores de la función en los x conocidos
## ts: valores en donde debe encontrarse la función interpolada
## retorna fs: valores de la función en los xs dados
function fs=interpole(t,f,ts)
  assert(size(t)==size(f));

  ts=ts(:); ## Asegúrese de que es un vector columna
  
  ## Encuentre las segundas derivadas
  fpp=findDerivs(t,f);

  ## ##################
  ## ## Problema 2.5 ##
  ## ##################

  ## >>> Ponga su solución aquí <<<
  
  ## Sugerencia: Puede serle muy útil el uso de la función 'lookup'
  ##             para encontrar cuál subintervalo utilizar.

  fs=zeros(size(ts));

	for iValue=1:length(ts)
		value = ts(iValue);

		i = lookup(t(2:length(t)),value)+2; # +1 para que primer indice sea 1
																				# +1 para que primer polinomio a utilizar sea el segundo
		if(i==length(t)+1) i=1; ## Condición xi+1 = x1
		endif
		
		iPrev = 1;
		if(i==1) iPrev=length(t); ## Punto anterior al inicial es el ultimo punto
		else iPrev=i-1;
		endif
		## Calculo de f(x)
		pt1 = fpp(iPrev)*( ((value-t(i))^3)*(1/(6*(t(iPrev)-t(i)))) );
		pt2 = fpp(i)*( ((value-t(iPrev))^3)*(1/(6*(t(i)-t(iPrev)))) );
		pt3 = ( (f(iPrev)/(t(iPrev)-t(i))) - ((fpp(iPrev)*(t(iPrev)-t(i)))/6) )*(value-t(i));
		pt4 = ( (f(i)/(t(i)-t(iPrev))) - ((fpp(i)*(t(i)-t(iPrev)))/6) )*(value-t(iPrev));
		fs(iValue) = pt1+pt2+pt3+pt4;
	endfor
	 
endfunction

## Depuración
figure(2,"name","Interpolación simple cerrada (depuración)");
x=[0,1,2,3];
f=[1,2,1,0.5];

hold off;
plot(x,f,'rx-;original;',"linewidth",2);

step=0.1;
xs=0:step:4-step;

fs=interpole(x,f,xs);

hold on;
plot(xs,fs,'bo-;interpolado;');

grid on;
xlabel("t");
ylabel("f(t)");

## El caso completo
N=10;
D = createData(N);
figure(1,"name","Interpolación 2D cerrada");
hold off;
plot(D(:,1),D(:,2),'rx-',"linewidth",2);

step=0.1;
t=0:step:N-step;
xs=interpole([0:N-1]',D(:,1),t);
ys=interpole([0:N-1]',D(:,2),t);

## ##################
## ## Problema 2.6 ##
## ##################

## >>> Ponga su solución aquí <<<



hold on;
plot(xs,ys,'bo-');
xlabel("x");
ylabel("y");
grid;
