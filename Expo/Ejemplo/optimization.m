1;          


# Area of trapezoid 
#		                      	_
#		\                   /		|
#		 \                 /		|
#			\ l           l /			| h
#		   \             /			|
#		    \___________/a    	|
#							b		
#
# With the constraint:
# 		b + 2*l = 7 [meters]
function z = Area(l,a)
	z = (7 - (l.*2) + (l.*cos(a))) * (l.*sin(a));
endfunction


# 
# Finds the maximun of a single variable function
# 
#
#
function x = goldenSection(l, a, xl, xu)
	epsilon = 10^-7;
	R = (sqrt(5)-1)/2;		# golden ratio
	d = R*(xu - xl);
	x1 = xl + d;
	x2 = xu - d;
	err = 1;
	xopt = 0;
	
	if(l==" ")
		f =@(x,y) Area(x,y);
		y = a;
	else
		f = @(x,y) Area(y,x);
		y = l;
	endif;
	
	while(err>epsilon)
		fx1 = f(x1,y);
		fx2 = f(x2,y);
		
		if(fx1 < fx2)
			xu = x1;
			x1 = x2;
			x2 = xu - R*(xu-xl);
			xopt = x2;
		else
			xl = x2;
			x2 = x1;
			x1 = xl + R*(xu-xl);
			xopt = x1;
		endif;
		
		err = (1-R)*abs((xu-xl)/xopt);
	end
	
	x = xopt;
	
endfunction;


#
# Finds the lenght and angle of a gutter
# that maximizes it's area, with the constraint
# that b + 2l = 6 meters
#
# Due to the constraint, the lenght of l
# can only be in the interval 0 < l < 3
#
function [l,a] = optimize()

	epsilon = 10^-7;
	
	current_a = 0.5;		# initial guess for the angle
	next_a = 0;
	current_l = goldenSection(" ", current_a, 0, 3.5, 1.5);
	next_l = 0;
	
	MAX_ITERATIONS = 100;
	
	aa=[];
	ll=[];
	
	for i = 1:MAX_ITERATIONS
		next_a = goldenSection(current_l, " ", 0, pi/2, pi/4);
		next_l = goldenSection(" ", next_a, 0, 3.5, 1.5);
		
		if(abs(next_a-current_a)<epsilon && abs(next_l-current_l)<epsilon)
			break;
		endif;
		
		current_a = next_a;
		current_l = next_l;
		aa(i)=next_a;
		ll(i)=next_l;	
	endfor;
	
	
	l = next_l;
	a = next_a;
	
	printf("Lenght of b: %f \n", 6-2*l);
	printf("Lenght of l: %f \n", l);
	printf("Angle: %f \n", a);
	printf("Area: %f \n", Area(l,a));

endfunction


function t = AreaPlot()

	l = linspace(0,3,100);
	a = linspace(0,pi/2,100);

	[ll,aa] = meshgrid(l,a);

	meshc(ll,aa, (ll.*sin(aa))*(6-(ll.*2)+(ll.*cos(aa))) );

endfunction






