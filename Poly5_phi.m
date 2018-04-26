function s_C1 = Poly5_phi(N,x1,x2,s1,s2,ds1,ds2,dds1,dds2,plotting)
%bx=linspace(x1,x2,N);
bx=[1:N]; 
dx=deriv(bx',1);
aa=0;
bb=0;
cc=0;
syms dd ee ff gg hh ii
%%x
eqn1 = aa*bx(1)^8+bb*bx(1)^7+cc*bx(1)^6+dd*bx(1)^5+ee*bx(1)^4+ff*bx(1)^3+gg*bx(1)^2+hh*bx(1)^1+ii*bx(1)^0 == s1;
eqn2 = aa*bx(end)^8+bb*bx(end)^7+cc*bx(end)^6+dd*bx(end)^5+ee*bx(end)^4+ff*bx(end)^3+gg*bx(end)^2+hh*bx(end)^1+ii*bx(end)^0 == s2;
%%spacing
eqn3 = 8*aa*bx(1)^7+7*bb*bx(1)^6+6*cc*bx(1)^5+5*dd*bx(1)^4+4*ee*bx(1)^3+3*ff*bx(1)^2+2*gg*bx(1)^1+1*hh*bx(1)^0+0*ii*bx(1)^0 == ds1/dx(1);
eqn4 = 8*aa*bx(end)^7+7*bb*bx(end)^6+6*cc*bx(end)^5+5*dd*bx(end)^4+4*ee*bx(end)^3+3*ff*bx(end)^2+2*gg*bx(end)^1+1*hh*bx(end)^0+0*ii*bx(end)^0 == ds2/dx(end);
%%second derivative of x
eqn6 = 7*8*aa*bx(1)^6+6*7*bb*bx(1)^5+5*6*cc*bx(1)^4+4*5*dd*bx(1)^3+3*4*ee*bx(1)^2+2*3*ff*bx(1)^1+1*2*gg*bx(1)^0+0*1*hh*bx(1)^0+0*ii*bx(1)^0 == dds1/dx(1);
eqn5 = 7*8*aa*bx(end)^6+6*7*bb*bx(end)^5+5*6*cc*bx(end)^4+4*5*dd*bx(end)^3+3*4*ee*bx(end)^2+2*3*ff*bx(end)^1+1*2*gg*bx(end)^0+0*1*hh*bx(end)^0+0*ii*bx(end)^0 == dds2/dx(end);

[As,Bs] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [dd,ee,ff,gg,hh,ii]);
Sol1 = linsolve(As,Bs);
s_C1=double(Sol1(1,1))*bx.^5+double(Sol1(2,1))*bx.^4+double(Sol1(3,1))*bx.^3+double(Sol1(4,1))*bx.^2+double(Sol1(5,1))*bx.^1+double(Sol1(6,1))*bx.^0;

if plotting=='t'
    figure
    plot(s_C1)
end
end

