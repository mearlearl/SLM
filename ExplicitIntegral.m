
function y = ExplicitIntegral(ff,LB,UB,dx)

    x = [LB:dx:UB];
    f = zeros(1,length(x));
    
    for i = 1:length(x)
    f(i) = ff(x(i));
    end
    
    m = length(f);
    y = 0;
    for i = 1:2:m-2
        y = y + dx*(f(i)+4*f(i+1)+f(i+2))/3;
    end    
    
end



