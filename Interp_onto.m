function Nc = Interp_onto(N,C,os,ns,plotting)
x_oc=C(:,1);
y_oc=C(:,3);

% Cubic Spline interpolation
x_nc=spline(os,x_oc,ns);
y_nc=spline(os,y_oc,ns);
% x_nc=pchip(os,x_oc,ns);
% y_nc=pchip(os,y_oc,ns);

x_nc(N)=x_oc(end);
y_nc(N)=y_oc(end);
x_nc(1)=x_oc(1);
y_nc(1)=y_oc(1);

Nc(:,1)=x_nc;
Nc(:,2)=y_nc;

if plotting=='t'
    figure
    plot(x_nc,y_nc,'-r')
    hold on
    plot(x_oc,y_oc,'-b')
end

end