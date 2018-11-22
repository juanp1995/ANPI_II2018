#***    Funcion que calcula la diferencia hacia atras           ***#
#***    para ello recibe el valor de h y x (1 en este ejemplo)  ***#
#***    y con las derivadas obtiene el valor.                   ***#

function [result] = diferenciaAtras(abc, X, Y, h)
funcion=(Y - (abc(1)*(X^2) + abc(2)*X + abc(3)))^2;
pderivadaA=2*(Y - (abc(1)*(X^2) + abc(2)*X + abc(3)))*(-X^2);
pderivadaB=2*(Y - (abc(1)*(X^2) + abc(2)*X + abc(3)))*(-X);
pderivadaC=2*(Y - (abc(1)*(X^2) + abc(2)*X + abc(3)))*(-1);
sderivadaA=2*(X^4);
sderivadaB=2*(X^2);
sderivadaC=2;


tempResultA=(funcion)-(pderivadaA*h)+((sderivadaA*(potencia(h,2)))/factorial(2));
tempResultB=(funcion)-(pderivadaB*h)+((sderivadaB*(potencia(h,2)))/factorial(2));
tempResultC=(funcion)-(pderivadaC*h)+((sderivadaC*(potencia(h,2)))/factorial(2));
result=[tempResultA,tempResultB,tempResultC];

endfunction