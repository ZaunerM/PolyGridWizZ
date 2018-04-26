function s = calc_s(x,y)
% dx=deriv(x,1);
% dy=deriv(y,1);
s(1)=0;
for i=2:size(x,1)
   dx = x(i,1)-x(i-1,1);
   dy = y(i,1)-y(i-1,1);
   ds = sqrt(dx^2+dy^2);
   s(i) = s(i-1)+ds;
end
end

