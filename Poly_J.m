function y_TE = Poly_J(N1,N2,sYs,sYc,sYe,dYs,dYc,dYe,ddYs,ddYc,ddYe,plotting)
%bx=[1:N1];
bx=linspace(sYs,sYc,N1);
dx=deriv(bx',1);
%Section 1 - Close to wall

syms dd ee ff gg hh ii
%x
eqn1 = dd*bx(1)^5+ee*bx(1)^4+ff*bx(1)^3+gg*bx(1)^2+hh*bx(1)^1+ii*bx(1)^0 == sYs;
eqn2 = dd*bx(end)^5+ee*bx(end)^4+ff*bx(end)^3+gg*bx(end)^2+hh*bx(end)^1+ii*bx(end)^0 == sYc; %0.7 %sYe;
%spacing
eqn3 = 5*dd*bx(1)^4+4*ee*bx(1)^3+3*ff*bx(1)^2+2*gg*bx(1)^1+1*hh*bx(1)^0+0*ii*bx(1)^0 == dYs/dx(1);
eqn4 = 5*dd*bx(end)^4+4*ee*bx(end)^3+3*ff*bx(end)^2+2*gg*bx(end)^1+1*hh*bx(end)^0+0*ii*bx(end)^0 == dYc/dx(end); %0.0032 %dYe;
%second derivative of x
eqn5 = 4*5*dd*bx(1)^3+3*4*ee*bx(1)^2+2*3*ff*bx(1)^1+1*2*gg*bx(1)^0+0*1*hh*bx(1)^0+0*ii*bx(1)^0 == ddYs/(dx(1)^2); %0
eqn6 = 4*5*dd*bx(end)^3+3*4*ee*bx(end)^2+2*3*ff*bx(end)^1+1*2*gg*bx(end)^0+0*1*hh*bx(end)^0+0*ii*bx(end)^0 == ddYc/(dx(end)^2); %1.2e-5 %ddYe; %6.5e-5
%3rd deriv
[As,Bs] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [dd,ee,ff,gg,hh,ii]);
Sol1 = linsolve(As,Bs);
y_TE1=double(Sol1(1,1)).*bx.^5+double(Sol1(2,1)).*bx.^4+double(Sol1(3,1)).*bx.^3+double(Sol1(4,1)).*bx.^2+double(Sol1(5,1)).*bx.^1+double(Sol1(6,1)).*bx.^0;

%%section2
%bx=[1:N2];
bx=linspace(sYc,sYe,N2);
dx=deriv(bx',1);
syms dd ee ff gg hh ii
%x
eqn1 = dd*bx(1)^5+ee*bx(1)^4+ff*bx(1)^3+gg*bx(1)^2+hh*bx(1)^1+ii*bx(1)^0 == sYc;
eqn2 = dd*bx(end)^5+ee*bx(end)^4+ff*bx(end)^3+gg*bx(end)^2+hh*bx(end)^1+ii*bx(end)^0 == sYe;
%spacing
eqn3 = 5*dd*bx(1)^4+4*ee*bx(1)^3+3*ff*bx(1)^2+2*gg*bx(1)^1+1*hh*bx(1)^0+0*ii*bx(1)^0 == dYc/dx(1);
eqn4 = 5*dd*bx(end)^4+4*ee*bx(end)^3+3*ff*bx(end)^2+2*gg*bx(end)^1+1*hh*bx(end)^0+0*ii*bx(end)^0 == dYe/dx(end);
%second derivative of x
eqn5 = 4*5*dd*bx(1)^3+3*4*ee*bx(1)^2+2*3*ff*bx(1)^1+1*2*gg*bx(1)^0+0*1*hh*bx(1)^0+0*ii*bx(1)^0 == ddYc/(dx(1)^2);
eqn6 = 4*5*dd*bx(end)^3+3*4*ee*bx(end)^2+2*3*ff*bx(end)^1+1*2*gg*bx(end)^0+0*1*hh*bx(end)^0+0*ii*bx(end)^0 == ddYe/(dx(end)^2); %6.5e-5
%3rd deriv
[As,Bs] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [dd,ee,ff,gg,hh,ii]);
Sol1 = linsolve(As,Bs);
y_TE2=double(Sol1(1,1)).*bx.^5+double(Sol1(2,1)).*bx.^4+double(Sol1(3,1)).*bx.^3+double(Sol1(4,1)).*bx.^2+double(Sol1(5,1)).*bx.^1+double(Sol1(6,1)).*bx.^0;

y_TE(1:size(y_TE1,2))=y_TE1;
y_TE(size(y_TE1,2)+1:size(y_TE1,2)+size(y_TE2,2)-1)=y_TE2(2:end);



if plotting=='t'
    figure
    plot(y_TE,deriv(y_TE',1))
    plot(deriv(y_TE',1))
%plot(deriv(deriv(y_TE',1),1))
end
end

