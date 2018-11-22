#***    Funcion que calcula la potencia de un numero   ***#
#***    para ello recibe el numero y la potencia a     ***#
#***    la cual se desea elevar.                       ***#

function [result] = potencia(num, exp)
cont=0;
result=1;
if(exp>=0)
while(cont~=exp)
result=result.*num;
cont+=1;
endwhile
else
while(cont~=exp)
result=result./num;
cont-=1;
endwhile
endif
endfunction