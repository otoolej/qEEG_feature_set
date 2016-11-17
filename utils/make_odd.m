function x=make_odd(x)
x=floor(x);
if(~rem(x,2)), x=x-1; end
