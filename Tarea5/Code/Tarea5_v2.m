1;

#
# Ordinary differential equation
#       y'=x*y^2
#
function u = f1(x,y)
	u = x*(y^2);
endfunction

#
# Stiff differential equation
#      y'=100-y
#
function u = f2(x,y)
	u = 100-y;
endfunction


#usage [x,y] = rungekutta4(f,xi,xf,y0,h)
#
# Solve ordinary differential equations using
# the Runge-Kutta method of 4th order
#
# f is the differential equation to be solve
#
# xi and xf are the initial and final points
#
# y0 initial value
#
# h is the step size
function [x,y] = rungekutta4(f,xi,xf,y0,h)


	current_x = xi;	
	current_y = y0;
	
	x_vector = zeros(1,1/h);
	y_vector = zeros(1,1/h);
	h_vector = zeros(1,1/h);
  i=1;
	
	while (current_x < xf) 

		#New ki values of formula.
		k1 = f(current_x,current_y);	
		k2 = f(current_x+(1/2)*h,current_y+(1/2)*k1*h);
		k3 = f(current_x+(1/2)*h,current_y+(1/2)*k2*h);
		k4 = f(current_x+h, current_y+k3*h);

		#Runge Kutta 4th degree iteration formula.
		next_y = current_y+(1/6)*(k1+2*k2+2*k3+k4)*h; 
		current_y = next_y;
		
		#Changes to the next x by adding the step to the current x.
		next_x = current_x+h;	
		current_x = next_x;

		x_vector(i)= current_x;
		y_vector(i)= current_y;
		h_vector(i) = h;
		i=i+1;
		
	end
	plot(x_vector,y_vector, 'linewidth', 1.5);
	
	#Plot configuration
  xlabel("x_{i}", "fontweight", "bold");
  ylabel("y(x_{i})", "fontweight", "bold");
  s = sprintf("Function evaluated in the interval %u \\leq x \\leq %u", xi, xf);
  title(s);
  box on;

	x=x_vector;				#Points where the function was evaluated
	y=y_vector; 	#Solution in the final point xf
endfunction



#usage [hi,yi] = rungekuttaError()
#
# Plots the error in the point x=1 in function
# of the step size h, of the ODE f1 = x*(y^2)
#
# Return values
#	hi Vector with the steps sizes
# yi Vector with the values of y(1)
#
function [hi,yi] = rungekuttaError()

	xi=0;							#Initial point
	xf=1;							#Final point
	y0=1;							#Initial condition
	h = 1/8;					#Initital step size
	realValue = 2;		#Real value of the function in x=1
	
	hi = zeros(1,8);
	yi = zeros(1,8);
	i=1;
	
	while(h>=(1/1024))
		
		[x,y] = rungekutta4(@f1,xi,xf,y0,h);
		hi(i) = h;
		yi(i) = abs(y(length(y))-realValue);
		h=h/2;
		i=i+1;
	end
	
	semilogy(hi,yi, 'linewidth', 1.5);
	
	#Plot configuration
	axis([min(hi), max(hi), min(yi), max(yi)]);
	#set(gca, 'XTick', [1/1024:127/5120:1/8]); # 5 ticks in the x-axis
  xlabel("h_{i}", "fontweight", "bold");
  ylabel("y(1)", "fontweight", "bold");
  s = sprintf("Value of the function y'=xy^2 evaluated in x=1 \n \
  \using different step sizes %s \\leq h \\leq %s","1/1024", "1/8");
  title(s);
  box on;

endfunction


#
#usage [t1,t2,t3]=stiffEquation()
#
# Plots the solution of the stiff differential
# equation y'=100-y using three different functions
#
# Shows the time required by each function and the
# number of points used.
#
function [t1,t2,t3] = stiffEquation()

	xi = 0;
	xf = 200;
	y0 = 5;
	h = xf/1000;
	
	rk_t1 = tic;
	[x,y] = rungekutta4(@f2, xi, xf, y0, h);
	rk_t2 = toc(rk_t1);
	
	ode45_t1 = tic;
	[x45,y45] = ode45(@f2, [xi, xf], y0);
	ode45_t2 = toc(ode45_t1);
	
	ode23_t1 = tic;
	[x23,y23] = ode23(@f2, [xi, xf], y0);
	ode23_t2 = toc(ode23_t1);
	
	subplot(2,1,1);
	hold on;
	plot(x,y, 'r', 'linewidth', 1.5);
	plot(x45,y45, 'g', 'linewidth', 1.5);
	plot(x23,y23, 'b', 'linewidth', 1.5);
	axis([0,200, 0,101]);
	
	xlabel("x_{i}", "fontweight", "bold");
  ylabel("y(x_{i})", "fontweight", "bold");
  s = sprintf("Solution of y'=100-y in the interval %u \\leq x \\leq %u", xi, xf);
  title(s);
  rk_leg = sprintf("Runge-Kutta\n Time: %u s \nPoints: %u", rk_t2, length(x));
  ode45_leg = sprintf("ode45\n Time: %u s \nPoints: %u", ode45_t2, length(x45));
  ode23_leg = sprintf("ode23\n Time: %u s \nPoints: %u", ode23_t2, length(x23));
  legend(rk_leg, ode45_leg, ode23_leg);
  legend('boxoff');
  legend('location', 'northeastoutside');
	hold off
	
	subplot(2,1,2);
	hold on;
	plot(x,y, 'r', 'linewidth', 1.5);			#runge-kutta plot
	plot(x45,y45, 'g', 'linewidth', 1.5);	#ode45 plot
	plot(x23,y23, 'b', 'linewidth', 1.5);	#ode23 plot
	
	axis([100, 200, 99.8, 100.2]);
	xlabel("x_{i}", "fontweight", "bold");
  ylabel("y(x_{i})", "fontweight", "bold");
  s = sprintf("Plot in the interval \n \
  \ x \\in [100,200] and y \\in [99.8, 100.2]");
  title(s);
  legend("rungekutta4", "ode45", "ode23");
  legend('boxoff');
  legend('location', 'northeastoutside');

endfunction




