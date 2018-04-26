function s_C1 = Poly5(s_tot,Ncl,NLE,Ncu,s1,scl,sLE,scu,s2,ds1,dsLE,ds2,plotting)
bx=s_tot;
aa=0;
bb=0;
cc=0;
syms aa bb cc dd ee ff gg hh ii
%%x
eqn1 = aa*bx(1)^8+bb*bx(1)^7+cc*bx(1)^6+dd*bx(1)^5+ee*bx(1)^4+ff*bx(1)^3+gg*bx(1)^2+hh*bx(1)^1+ii*bx(1)^0 == s1;
eqn2 = aa*bx(Ncu)^8+bb*bx(Ncu)^7+cc*bx(Ncu)^6+dd*bx(Ncu)^5+ee*bx(Ncu)^4+ff*bx(Ncu)^3+gg*bx(Ncu)^2+hh*bx(Ncu)^1+ii*bx(Ncu)^0 == scu;
eqn3 = aa*bx(NLE)^8+bb*bx(NLE)^7+cc*bx(NLE)^6+dd*bx(NLE)^5+ee*bx(NLE)^4+ff*bx(NLE)^3+gg*bx(NLE)^2+hh*bx(NLE)^1+ii*bx(NLE)^0 == sLE;
eqn4 = aa*bx(Ncl)^8+bb*bx(Ncl)^7+cc*bx(Ncl)^6+dd*bx(Ncl)^5+ee*bx(Ncl)^4+ff*bx(Ncl)^3+gg*bx(Ncl)^2+hh*bx(Ncl)^1+ii*bx(Ncl)^0 == scl;
eqn5 = aa*bx(end)^8+bb*bx(end)^7+cc*bx(end)^6+dd*bx(end)^5+ee*bx(end)^4+ff*bx(end)^3+gg*bx(end)^2+hh*bx(end)^1+ii*bx(end)^0 == s2;
%%spacing
eqn6 = 8*aa*bx(1)^7+7*bb*bx(1)^6+6*cc*bx(1)^5+5*dd*bx(1)^4+4*ee*bx(1)^3+3*ff*bx(1)^2+2*gg*bx(1)^1+1*hh*bx(1)^0+0*ii*bx(1)^0 == ds1;
eqn7 = 8*aa*bx(NLE)^7+7*bb*bx(NLE)^6+6*cc*bx(NLE)^5+5*dd*bx(NLE)^4+4*ee*bx(NLE)^3+3*ff*bx(NLE)^2+2*gg*bx(NLE)^1+1*hh*bx(NLE)^0+0*ii*bx(NLE)^0 == dsLE;
eqn8 = 8*aa*bx(end)^7+7*bb*bx(end)^6+6*cc*bx(end)^5+5*dd*bx(end)^4+4*ee*bx(end)^3+3*ff*bx(end)^2+2*gg*bx(end)^1+1*hh*bx(end)^0+0*ii*bx(end)^0 == ds2;

[As,Bs] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8], [aa,bb,cc,dd,ee,ff,gg,hh,ii]);
Sol1 = linsolve(As,Bs);
s_C1=double(double(Sol1(1,1))*bx.^8+double(Sol1(2,1))*bx.^7+double(Sol1(3,1))*bx.^6+double(Sol1(4,1))*bx.^5+double(Sol1(5,1))*bx.^4+double(Sol1(6,1))*bx.^3+double(Sol1(7,1))*bx.^2+double(Sol1(8,1))*bx.^1+double(Sol1(9,1)));

if plotting=='t'
    figure
    plot(s_C1)
end
end

