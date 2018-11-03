function x = interpolation(l, a, x1, x2, x3)

	epsilon = 10^-15;
	x4 = x3+1;
	y=0;
	
	if(l==" ")
		f =@(x,y) Area(x,y);
		y = a;
	else
		f = @(x,y) Area(y,x);
		y = l;
	endif;
	
	while(abs(x4-x3)>epsilon)
		num = f(x1,y)*(x2^2 - x3^2) + f(x2,y)*(x3^2 - x1^2) + f(x2,y)*(x1^2 - x2^2);
		den = f(x1,y)*(x2 - x3) + f(x2,y)*(x3 - x1) + f(x3,y)*(x1 - x2);
		
		x4 = 0.5*(num/den);
		
		if(abs(x1-x4) < abs(x2-x4))
			x2 = x3;
			x3 = x4;
		else
			x1 = x3;
			x3 = x4;
		endif
	end
	
	x = x4;

endfunction