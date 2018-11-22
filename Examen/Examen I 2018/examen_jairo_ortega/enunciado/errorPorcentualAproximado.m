#***    Funcion que calcula el error porcentual aproximado   ***#
#***    para un metodo iterativo, para ello recibe la        ***#
#***    aproximacion actual y la aproximacion anterior.      ***#

function [ea] = errorPorcentualAproximado (aproxAct, aproxAnt)
ea=((aproxAct-aproxAnt)/(aproxAct))*100;
endfunction