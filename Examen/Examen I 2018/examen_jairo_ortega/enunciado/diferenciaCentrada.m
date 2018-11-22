#***    Funcion que calcula la diferencia centrada                  ***#
#***    para ello recibe el valor de h y x (1 en este ejemplo)      ***#
#***    y con las diferencias de atras y adelante obtiene el valor. ***#

function [result] = diferenciaCentrada(abc, X, Y, h)
difAd=diferenciaAdelante(abc, X, Y, h);
difAt=diferenciaAtras(abc, X, Y, h);
result=(difAd-difAt)/(2*h);
endfunction