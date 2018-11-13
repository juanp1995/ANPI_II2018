1;          


# Area of trapezoid 
# A = (b + B)/2 * h
# 
#	d = l*cos(a)
# h = l*sin(a)
# B = b + 2*d
#						B	
#	__--------------------__  _ 
#		\                   /		|
#		 \                 /		|
#			\ l           l /			| h
#		   \             /			|
#		    \___________/a    	|
#	--d----			b		
#
# With the constraint:
# 		b + 2*l = 7 [meters]
# 	--> b = 7 - 2*l




#
# Area of the gutter in terms
# of the lenght of the sides 
# and the angle on inclination
#
function z = Area(l,a,s)
	b = s - l.*2;
	z = (b + l.*cos(a))*l.*sin(a);
	#z = (7 - (l.*2) + (l.*cos(a))) * (l.*sin(a));

endfunction


# 
# Finds the maximun of a single variable function
# 
#
#
function x = goldenSection(l, a, s, xl, xu)
	epsilon = 10^-7;
	R = (sqrt(5)-1)/2;		# golden ratio
	d = R*(xu - xl);
	x1 = xl + d;
	x2 = xu - d;
	err = 1;
	xopt = 0;
	
	if(l==" ")
		f =@(x,y) Area(x,y, s);
		y = a;
	else
		f = @(x,y) Area(y,x, s);
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
function [l,a] = optimize(s)

	epsilon = 10^-7;
	cons_l = s/2;
	
	##current_a = 0.5;		# initial guess for the angle
	current_a = 10;
	next_a = 0;
	current_l = goldenSection(" ", current_a, s, 0, cons_l, cons_l/2);
	next_l = 0;
	
	MAX_ITERATIONS = 500;
	
	aa=[];
	ll=[];
	
	for i = 1:MAX_ITERATIONS
		next_a = goldenSection(current_l, " ", s, 0, pi/2, pi/4);
		#next_a = goldenSection(current_l, " ", 0, 90, 45);
		next_l = goldenSection(" ", next_a, s, 0, cons_l, cons_l/2);
		
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
	
	printf("Lenght of b: %f \n", s-2*l);
	printf("Lenght of l: %f \n", l);
	printf("Angle: %f \n", a*(180/pi));
	printf("Area: %f \n", Area(l,a,s));

endfunction


function t = AreaPlot()

	l = linspace(-10,10,100);
	a = linspace(-10,10,100);

	[ll,aa] = meshgrid(l,a);

	mesh(ll,aa, (ll.*sin(aa))*(6-(ll.*2)+(ll.*cos(aa))) );

endfunction






