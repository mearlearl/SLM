function xzero = newton(f,fp,x0,tol)
%NEWTON find roots of equation f
err = 1;
while(err>tol)
    xzero = x0 - f(x0)/fp(x0);
    err = abs(xzero-x0)/xzero;
    x0 = xzero;
end
end

