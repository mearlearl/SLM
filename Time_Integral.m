function y = Time_Integral(f,dx)
  
    sz = size(f);
    m = sz(2); 
    y = 0;
    for i = 1:2:m-2
        y = y + dx*(f(:,i)+4*f(:,i+1)+f(:,i+2))/3;
    end    
    
end

