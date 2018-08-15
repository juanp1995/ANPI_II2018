1;


# usage [v, ev, ea, n] = anpi_erf (x, c)
#
# Cálculo estimado de la función de error utilizando la serie de Taylor
#
# v   = Valor estimado
# ev  = Error relativo porcentual verdadero
# ea  = Error relativo porcentual aproximado
# n   = Número de terminos de la serie utilizados para el cálculo
#
# x   = Valor real de entrada
# c   = Número de cifras significativas, debe ser menor a 16 (precision doble)
#       (por defecto es 5)
#
#
function [v, ev, ea, n] = anpi_erf (x, c)
  
  if (nargin != 1),
    if(c > 15),
    error("Cifras significativas deben ser menor de 15 (precisión doble)");
    end;
  else
    c=5;
  end;  
  
  v=1;
  ev = n =0;
  ea=100;
  es = (0.5*10^(2-c));
  rv = erf(x); #Valor real
  output_precision(c+1);
  
  aprox = 1;
  aprox_ant = 0;
  v_ant = 0;
  const = (2/sqrt(pi));
  
  while(abs(ea) > es);
    aprox =  (((((-1)^n)*x^(2*n+1))/(factorial(n)*(2*n+1))) + aprox_ant);
    v = (const * aprox);
    v_ant = (const * aprox_ant);
    ea = (((v - v_ant)/v) * 100);
    n=n+1;
    aprox_ant = aprox;
  end;
  
  ev = ((1-(v/rv)) * 100);
 
endfunction


# usage anpi_quadratic(a, b, c) 
#
#Cálculo de raíces de ecuaciónes cuadráticas de la forma 
#         ax^2 + bx + c = 0
#
# Uso;  1) anpi_quadratic(a, b, c);
#           Imprime en consola los valores de las raíces calculados con
#           la función tradicional y alternativa, tanto en precisión
#           simple como doble.
#
#       2) [x1, x2, x1_alt, x2_alt] = anpi_quadratic(a, b, c);
#           No imprime los valores en consola, en los parametros de salida
#           se retornan los valores de las raíces calculados con la función
#           tradicional y alternativa en precisión doble.
#
# x1      = Primera raíz con la formula tradicional
# x2      = Segunda raíz con la formula tradicional
# x1_alt  = Primera raíz con la formula alternativa
# x2_alt  = Segunda raíz con la formula tradicional
#
# a       = Coeficiente cuadratico
# b       = Coeficiente lineal
# c       = Coeficiente independiente
#
#
function [x1, x2, x1_alt, x2_alt] = anpi_quadratic (a, b, c)
    
  #Cálculo
  x1 = double((-b + sqrt(b^2-4*a*c))/2*a);
  x2 = double((-b - sqrt(b^2-4*a*c))/2*a); 
  
  x1_alt = double(((-2*c)/(b+sqrt(b^2-4*a*c))));
  x2_alt = double(((-2*c)/(b-sqrt(b^2-4*a*c))));
       

  if (nargout == 0),
    #Salida
    printf("Resultado fórmula tradicional\n");
    printf("Precisión simple: \n");
    format long e;
    printf("X1 ="); disp(single(x1));
    printf("X2 ="); disp(single(x2));
    printf("Precisión doble: \n");
    printf("X1 ="); disp(x1);
    printf("X2 ="); disp(x2);
    printf("\n");
    
    printf("Resultado formula alternativa\n");
    printf("Precisión simple: \n");
    printf("X1 ="); disp(single(x1_alt));
    printf("X2 ="); disp(single(x2_alt));
    printf("Precisión doble: \n");
    printf("X1 ="); disp(x1_alt);
    printf("X2 ="); disp(x2_alt);
  end;  
  
endfunction