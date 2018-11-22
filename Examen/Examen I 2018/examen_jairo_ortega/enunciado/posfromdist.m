#!/usr/bin/octave-cli

## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## I Semestre 2018
## Examen Final

## PROBLEMA 2

## NOMBRE: Jairo Ortega Calderón
## CARNE: 2014043224


## Cargue los datos de distancias
##
## Cada columna de D tiene las distancias a un sensor fijo en el
## espacio. 
D=load("-ascii","dists.dat");

## Las posiciones de los emisores están dadas por las _columnas_
## de la matriz S
E=load("-ascii","emisors.dat");

## Función calcula posición a partir de distancias a emisores.
##
## dists:     Distancias del punto a cada emisor como vector fila.
## emisorPos: Posiciones de los emisores, cada emisor en una columna.
## option:    Si solo hay 3 emisores, cuál de las dos soluciones retorna
## La posición i de dists contiene la distancia al emisor en la
## columna i.
##
function p=calcPosition(dists,emisorPos,option=1)

  ## Número de emisores usados
  dim = min(columns(emisorPos),length(dists));

  ## ##################
  ## ## Problema 2.1 ##
  ## ##################
  ## Construya la matriz M
  M  = zeros(dim,4); ## CAMBIE ESTO!!
  contM=1;
  while(contM<=dim)
    M(contM,1) = 1;
    M(contM,2) = -2 * emisorPos(1,contM);
    M(contM,3) = -2 * emisorPos(2,contM);
    M(contM,4) = -2 * emisorPos(3,contM);
    contM=contM+1;
  endwhile
  ## ##################
  ## ## Problema 2.2 ##
  ## ##################
  ## Construya el vector b
  
  b = zeros(dim,1); 
  contB=1;
  while(contB<=dim)
  b(contB,1) = (dists(1,contB) ^ 2)-((emisorPos(1,contB))^2 + (emisorPos(2,contB))^2 + (emisorPos(3,contB))^2); 
  contB=contB+1;
  endwhile

  ## ##################
  ## ## Problema 2.3 ##
  ## ##################

  ## Calcule la matriz seudo-inversa utilizando SVD
  
  [U,S,V] = svd(M);
  Ur = resize(U,dim,4);
  Vr = resize(V,4,4);
  Sr = resize(S,4,4);
  iM =  Vr * inv(Sr) * Ur';
  
  ## Verifique que iM y pinv(M) son lo mismo
  if (norm(iM-pinv(M),"fro") > 1e-6)
    error("Matriz inversa calculada con SVD incorrecta");
  endif

  ## ##################
  ## ## Problema 2.4 ##
  ## ##################
  
  ## Calcule la solución particular
  hatp=zeros(4,1); ## CAMBIE ESTO!!
  hatp = iM * b;  
  
  ## El caso de 3 dimensiones tiene dos posibles soluciones:
  if (dim==3)

    ## ##################
    ## ## Problema 2.5 ##
    ## ##################

    ## Con 3 emisores, calcule las dos posibles posiciones

    ## PONGA SU CÓDIGO AQUÍ
    
    hatp2 = pinv(M) * b;
    vectNulo = [Vr(1,4);Vr(2,4);Vr(3,4);Vr(4,4)];
    a = vectNulo(2) ^ 2 + vectNulo(3) ^ 2 + vectNulo(4) ^ 2;
    b = (2 * vectNulo(2) * hatp2(2)+ 2 * vectNulo(3) * hatp2(3) + 2 * vectNulo(4) * hatp2(4) - vectNulo(1));
    c = ( hatp2(2) ^ 2 + hatp2(3) ^ 2 + hatp2(4) ^ 2 - hatp2(1));
    
    if(option==1)
     sol  = (-b + sqrt(b^2 - 4*a*c)) / (2 * a);
    else
      sol = (-b - sqrt(b^2 - 4*a*c)) / (2 * a);
    endif 
    p = [hatp2(1) + sol*vectNulo(1);hatp2(2) + sol*vectNulo(2); hatp2(3) + sol*vectNulo(3); hatp2(4) + sol*vectNulo(4)];
    
  else 
    ## ##################
    ## ## Problema 2.6 ##
    ## ##################

    ## Caso general
    
    ## PONGA SU CÓDIGO AQUÍ
    p = hatp;
  endif
  
endfunction

## Calcule la trayectoria de posiciones para una matriz de
## distancias que contiene en cada fila el vector de distancias
## a cada emisor.
## 
## Sea N el número de datos y S el número de emisores
## dists:     matriz de distancias de tamaño N x S
## emisorPos: matrix de tamaño 3xS
## option:    si S==3, entonces cuál de las dos soluciones devolver
##            option=1 implica usar + en cuadrática, y otra cosa usa -
function p=calcPositions(dists,emisorPos,option=1)
  n=min(columns(emisorPos),columns(dists));
  k=rows(dists);

  ## Reserve memoria para todos los puntos en la trayectoria
  p=zeros(k,4);
  
  ## Para todos los puntos en la trayectoria
  for i=1:k
    ## Calcule la posición del punto, dadas las distancias a los emisores
    p(i,:)=calcPosition(dists(i,1:n),emisorPos(:,1:n),option)';
  endfor

endfunction

## Función calcula distancia de un punto a todos los emisores
## pos:       Posición del objeto (en las filas)
## emisorPos: Posiciones de los emisores, en las columnas.
function d=calcDistances(pos,emisorPos)
  ## ##################
  ## ## Problema 2.7 ##
  ## ##################

  ## Calcule la distancia a cada emisor, dada la posición "pos" del
  ## objeto y las posiciones de los emisores.
  
  ## PONGA SU CÓDIGO AQUÍ
  
  for i = 1:rows(pos)
    for j = 1:5
      xe = emisorPos(1,j);
      ye = emisorPos(2,j);
      ze = emisorPos(3,j);
      d(i,j) = sqrt((pos(i,2)-xe)^2 + (pos(i,3)-ye)^2 + (pos(i,4)-ze)^2); # distancia
    endfor
  endfor
  d
endfunction


## #######################
## ## Pruebas unitarias ##
## #######################

## Calculando una posición a partir de 3 sensores
display("Usando 3 emisores\n");
p=calcPosition(D(1,1:3),E)

## Calculando una posición a partir de 4 sensores
display("Usando 4 emisores\n");
p=calcPosition(D(1,1:4),E)

## Calculando una posición a partir de 5 sensores
display("Usando 5 emisores\n");
p=calcPosition(D(1,1:5),E)

## ########################
## ## Número de emisores ##
## ########################
n=5; ## Cambie aquí el número de emisores a usar

## Calcule las posiciones a partir de las distancias
p=calcPositions(D(:,1:n),E);

## Distancias a posiciones estimadas
ed=calcDistances(p,E);

## Errores
err=(ed-D);

## Promedio de errores
sum(err)/rows(err)

## Muestre los errores de las distancias a los primeros tres emisores
figure(1,"name","Errors")
hold off;
plot(err(:,1),"r");
hold on;
plot(err(:,2),"b");
plot(err(:,3),"k");

## Grafique las posiciones encontradas
figure(2,"name","Trayectoria estimada y emisores");
hold off;

## Grafique la trayectoria en azul
plot3(p(:,1),p(:,2),p(:,3),'b');
xlabel("x");
ylabel("y");
zlabel("z");

hold on;

if (n==3)  
  ## Si n==3 hay otros posibles valores para la trayectoria
  p=calcPositions(D(:,1:n),E,2);
  plot3(p(:,1),p(:,2),p(:,3),'r');
  
  ## A la fuerza use n=5 para ver la trayectoria verdadera
  p=calcPositions(D(:,1:5),E);
  plot3(p(:,1),p(:,2),p(:,3),'go');
endif

## Grafique los emisores en postes
plot3([E(1,1:n);E(1,1:n)],[E(2,1:n);E(2,1:n)],[E(3,1:n);zeros(1,n)],"linewidth",3);
plot3(E(1,1:n),E(2,1:n),E(3,1:n),'rx');
