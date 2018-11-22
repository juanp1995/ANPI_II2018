#***    Funcion que calcula el factorial de un numero   ***#

function [result]= factorial(n)
result=1;
while(n~=0)
result*=n;
n=n-1;
endwhile
endfunction

