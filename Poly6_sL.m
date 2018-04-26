function s_C1 = Poly6_sL(N,bx,y1,y2,ds1,ds2,dds1,dds2,ddds2,plotting)

syms dd ee ff gg hh ii
%%x
eqn1 =dd*  bx(1)^5+ee*  bx(1)^4+ff*  bx(1)^3+gg*  bx(1)^2+hh*  bx(1)^1+ii == y1;
eqn2 =dd*bx(end)^5+ee*bx(end)^4+ff*bx(end)^3+gg*bx(end)^2+hh*bx(end)^1+ii == y2;
%%spacing
eqn3 = 5*dd*  bx(1)^4+4*ee*  bx(1)^3+3*ff*  bx(1)^2+2*gg*  bx(1)^1+1*hh == ds1;
eqn4 = 5*dd*bx(end)^4+4*ee*bx(end)^3+3*ff*bx(end)^2+2*gg*bx(end)^1+1*hh == ds2;
%%second derivative of x
eqn5 = 4*5*dd*  bx(1)^3+3*4*ee*  bx(1)^2+2*3*ff*  bx(1)^1+1*2*gg == dds1;
eqn6 = 4*5*dd*bx(end)^3+3*4*ee*bx(end)^2+2*3*ff*bx(end)^1+1*2*gg == dds2;
%%3rd deriv
%eqn6 = 4*5*6*cc*bx(1)^3+3*4*5*dd*bx(1)^2+2*3*4*ee*bx(1)^1+1*3*ff == ddds2;

[As,Bs] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [dd,ee,ff,gg,hh,ii]);
Sol1 = linsolve(As,Bs); %mldivide(As,Bs);
%Sol1  = vpa(Sol1_);
s_C1=double(Sol1(1,1))*bx.^5+double(Sol1(2,1))*bx.^4+double(Sol1(3,1))*bx.^3+double(Sol1(4,1))*bx.^2+double(Sol1(5,1))*bx.^1+double(Sol1(6,1));

if plotting=='t'
    figure
    plot(bx,s_C1)
end
end

