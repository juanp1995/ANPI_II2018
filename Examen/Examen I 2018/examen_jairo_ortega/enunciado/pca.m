## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## I Semestre 2018
## Examen Final

## PROBLEMA 3

## NOMBRE: Jairo Ortega Calderón
## CARNE: 2014043224

3;

## Cargue los datos 
X=load("-ascii","pcadata.dat");
N=columns(X);

################################################
## Problema 3.1                               ##
## Grafique los datos                         ##
################################################
figure(1);
plot3(X(1,1:end),X(2,1:end),X(3,1:end),'x');
grid on;
title('Matriz X')
xlabel("X")
ylabel("Y")
zlabel("Z")
hold off;

################################################
## Problema 3.2                               ##
## Calcule la media de los datos              ##
################################################

# El valor medio se calcula de la siguiente manera:
# u = (1/n)sum_{i=1}^{n}{Xi}
function [result] = valorMedio(N, X)
  totX=0;
  totY=0;
  totZ=0;
  cont=1;
  while(cont<=N)
  totX= totX + X(1,cont);
  totY= totY + X(2,cont);
  totZ =totZ + X(3,cont);
  cont+=1;
  endwhile
  Ux=(1/N)*(totX);
  Uy=(1/N)*(totY);
  Uz=(1/N)*(totZ);
  result = [Ux,Uy,Uz];
endfunction

################################################
## Problema 3.3                               ##
## Muestre la media en rojo                   ##
################################################

[U]=valorMedio(columns(X), X);
figure(2);
plot3(X(1,1:end),X(2,1:end),X(3,1:end),'x');
hold on;
plot3(U(1),U(2),U(3), 'or');
hold on;
grid on;
title('Matriz X con punto medio');
xlabel("X")
ylabel("Y")
zlabel("Z")
hold off;


################################################
## Problema 3.4                               ##
## Calcule los datos sin media                ##
################################################
[U]=valorMedio(columns(X), X);
xMed=[];
yMed=[];
zMed=[];
cont=1;
while(cont<=N)
xMed=[xMed;(X(1,cont)-U(1))];
yMed=[yMed;(X(2,cont)-U(2))];
zMed=[zMed;(X(3,cont)-U(3))];
cont+=1;
endwhile
XM0=[xMed, yMed, zMed]';
figure(3);
plot3(XM0(1,1:end),XM0(2,1:end),XM0(3,1:end),'x');
hold on;
grid on;
title('Matriz X de Media Cero');
xlabel("X")
ylabel("Y")
zlabel("Z")
hold off;


################################################
## Problema 3.5                               ##
## Calcule la matriz de covarianza            ##
################################################
# La matriz de covarianza se calcula de la siguiente manera:
# MatCov = (1/n)*X*X^{T}

XMatCov=(XM0*XM0').*(1/N)

################################################
## Problema 3.6                               ##
## Encuentre los eigenvalores y eigenvectores ##
################################################

[eigenvectores, eigenvalores]=eig(XMatCov)

################################################
## Problema 3.7                               ##
## Reordene los eigenvectores para PCA        ##
################################################

eigVecComponentesPrincipales=eigenvectores'

########################################################################
## Problema 3.8                                                       ##
## Cuáles son los ejes principales y qué varianza tiene los datos     ##
########################################################################

EjePrincipal1=eigVecComponentesPrincipales(1, 1:end)
EjePrincipal2=eigVecComponentesPrincipales(2, 1:end)
EjePrincipal3=eigVecComponentesPrincipales(3, 1:end)

varianzaEje1=eigenvalores(1,1)
varianzaEje2=eigenvalores(2,2)
varianzaEje3=eigenvalores(3,3)

########################################################################
## Problema 3.9                                                       ##
## Calcule la proyección de los datos al plano engendrado por los dos ##
## eigenvectores                                                      ##
########################################################################

#Proyección: Y=PX
Y=[EjePrincipal1;EjePrincipal2]*XM0;
## Grafique la proyección
figure(4);
plot(Y(1,1:end),Y(2,1:end), '.');
hold on;
grid on;
title('Puntos proyectados sobre el plano');
xlabel("X")
ylabel("Y")
hold off;

########################################################################
## Problema 3.10                                                      ##
## Calcule los datos reconstrudos a partir de los datos proyectados   ##
########################################################################
[U]=valorMedio(columns(X), X);
Xaprox=([EjePrincipal1;EjePrincipal2]'*Y)+U';

figure(5);
plot3(X(1,1:end),X(2,1:end),X(3,1:end),'x');
hold on;
plot3(Xaprox(1,1:end),Xaprox(2,1:end),Xaprox(3,1:end), 'sm');
hold on;
grid on;
title('Matriz reconstruida');
xlabel("X")
ylabel("Y")
zlabel("Z")
hold off;
