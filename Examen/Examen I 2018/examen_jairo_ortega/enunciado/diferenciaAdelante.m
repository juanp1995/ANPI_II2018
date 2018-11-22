#***    Funcion que calcula la diferencia hacia adelante        ***#
#***    para ello recibe el valor de h y x (1 en este ejemplo)  ***#
#***    y con las derivadas obtiene el valor.                   ***#

function [result] = diferenciaAdelante(abc, x, y, h)
funcion = (y - (abc(1)*(x^2) + abc(2)*x + abc(3)))^2;
pderivadaA=2*(y - (abc(1)*(x^2) + abc(2)*x + abc(3)))*(-x^2);
pderivadaB=2*(y - (abc(1)*(x^2) + abc(2)*x + abc(3)))*(-x);
pderivadaC=2*(y - (abc(1)*(x^2) + abc(2)*x + abc(3)))*(-1);
sderivadaA=2*(x^4);
sderivadaB=2*(x^2);
sderivadaC=2;

tempResultA=(funcion)+(pderivadaA*h)+((sderivadaA*(potencia(h,2)))/factorial(2));
tempResultB=(funcion)+(pderivadaB*h)+((sderivadaB*(potencia(h,2)))/factorial(2));
tempResultC=(funcion)+(pderivadaC*h)+((sderivadaC*(potencia(h,2)))/factorial(2));
result=[tempResultA,tempResultB,tempResultC];
endfunction
