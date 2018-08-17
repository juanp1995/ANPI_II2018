1;
#usage [fd_v] = derivada(x)
#
# Cálculo de la derivada de la función
#       f(x) = sin(e^(x^2))
#
# utilizando los métodos de diferenciación
#   * hacia atrás 
#   * hacia adelante
#   * centrada
#
# Muestra el gráfico del error relativo fraccional verdadero
# en función del tamaño de paso
#
#
function [fd_v] = derivada (x)

  muestras_decada = 100; 
  lambda = 10^(1/(-muestras_decada)); #Factor de reducción para 100 muestras/decada
 
  total_muestras = log10(10^0/10^-15)*muestras_decada;
  h = 1;            #Tamaño de paso
  x_sig = x + h;
  x_ant = x - h;

  f = @(z) sin(e^(z^2));                  #Función a derivar
  fd_v = (2*x*(e^(x^2)))*(cos(e^(x^2)));  #Valor real de la derivada
  

  #Vectores con los errores relativos de cada diferenciación
  err_ad = zeros(1,total_muestras);
  err_at = zeros(1,total_muestras);
  err_ce = zeros(1,total_muestras);
  
  #Vector con los tamaños de paso
  hi = zeros(1, total_muestras);
  i = 1;
  
  while (h > 10^-15)
    # Diferenciación hacia adelante
    fd_aprox_ad = (f(x_sig) - f(x)) / h;
    err_aprox_ad = abs((fd_v - fd_aprox_ad) / fd_v);
    
    err_ad(i) = err_aprox_ad;
    
    ## Diferenciación hacia atrás
    fd_aprox_at = (f(x) - f(x_ant)) / h;
    err_aprox_at = abs((fd_v - fd_aprox_at) / fd_v );
   
    err_at(i) = err_aprox_at; 
    
    ## Diferenciación centrada
    fd_aprox_ce = (f(x_sig) - f(x_ant)) / (2*h);
    err_aprox_ce = abs((fd_v - fd_aprox_ce) / fd_v );
    
    err_ce(i) = err_aprox_ce;
    
    hi(i) = h;
    i = i+1;
    ## Siguiente paso
    h = lambda*h;
    x_sig = x+h;
    x_ant = x-h;
  end;  
  
  #Muestra de datos en el gráfico
  hold on;
  loglog(hi, err_ad, 'r', 'linewidth', 1.5);
  loglog(hi, err_at, 'g', 'linewidth', 1.5);
  loglog(hi, err_ce, 'b', 'linewidth', 1.5);
  
  #Configuración del gráfico
  xlabel("h_i", "fontweight", "bold");
  ylabel("E_{rel}", "fontweight", "bold");
  s = sprintf("Error relativo fraccional vs tamaño de paso h \n x=%u", x);
  title(s);
  legend("adelante", "atras", "centrada");
  legend('boxon');
  legend('location', 'northwest');
  box on;
  
endfunction
