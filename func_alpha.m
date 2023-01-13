function result = func_alpha(a0,a1,t)

result = a1 +   (a0-a1)*( (1-t)- sin(2*pi*(1-t))/(2*pi)   );