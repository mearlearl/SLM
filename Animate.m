for j = 1:2*N
   p=plot(x,u(:,j),'r-');
   p.LineWidth = 3;
   axis([0 1 min(min(u)) max(max(u))+1]);
   pause(0.05)
end    
pause(0.5);
surf(u)