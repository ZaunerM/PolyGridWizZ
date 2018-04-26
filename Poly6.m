    function s_C1 = Poly6(N,s1,s2,ds1,ds2,dds1,dds2,ddds2,plotting)
bx=[1:N]; 
%bx=linspace(x1,x2,N);
dx=1; %deriv(bx',1);

syms cc dd ee ff gg hh ii
%%x
eqn1 =       cc * bx(1)^6 +      dd * bx(1)^5 +      ee * bx(1)^4 +      ff * bx(1)^3 +    gg * bx(1)^2 +   hh * bx(1)^1 +ii == s1;
eqn2 =       cc*bx(end)^6 +      dd*bx(end)^5 +      ee*bx(end)^4 +      ff*bx(end)^3 +    gg*bx(end)^2 +   hh*bx(end)^1 +ii == s2;
%%spacing
eqn3 =     6*cc * bx(1)^5 +    5*dd * bx(1)^4 +    4*ee * bx(1)^3 +    3*ff * bx(1)^2 +  2*gg * bx(1)^1 +1*hh           +0  == ds1/dx(1);
eqn4 =     6*cc*bx(end)^5 +    5*dd*bx(end)^4 +    4*ee*bx(end)^3 +    3*ff*bx(end)^2 +  2*gg*bx(end)^1 +1*hh           +0  == ds2/dx(end);
%%second derivative of x
eqn7 =   5*6*cc * bx(1)^4 +  4*5*dd * bx(1)^3 +  3*4*ee * bx(1)^2 +  2*3*ff * bx(1)^1 +1*2*gg           +0              +0  == dds1/(dx(1)^2);
eqn5 =   5*6*cc*bx(end)^4 +  4*5*dd*bx(end)^3 +  3*4*ee*bx(end)^2 +  2*3*ff*bx(end)^1 +1*2*gg           +0              +0  == dds2/(dx(end)^2);
%%3rd deriv
eqn6 = 4*5*6*cc * bx(1)^3 +3*4*5*dd * bx(1)^2 +2*3*4*ee * bx(1)^1 +1*2*3*ff           +0                +0              +0  == ddds2/dx(1);

[As,Bs] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7], [cc,dd,ee,ff,gg,hh,ii]);

%Sol1 = linsolve(As,Bs);
Sol1 = mldivide(As,Bs);

s_C1_=double(Sol1(1,1))*bx.^6+double(Sol1(2,1))*bx.^5+double(Sol1(3,1))*bx.^4+double(Sol1(4,1))*bx.^3+double(Sol1(5,1))*bx.^2+double(Sol1(6,1))*bx.^1+double(Sol1(7,1))*bx.^0;
s_C1=s_C1_;
if plotting=='t'
    figure
    %plot(deriv(deriv(s_C1',1)./dx,1)./dx)
    hold on
     plot(deriv(s_C1_',1),'r')
%    plot(s_C1_','r')
end

% r1=linspace(0,(N-1)*ds1,N);
% r2=linspace(s2,s2-(N-1)*ds2,N);
% 
% lam=[1:Nb];
% a=-2/( 2*((1-Nb^3)-3*(1-Nb)) - ((1-Nb^2)-2*(1-Nb))*3*(1+Nb) );
% b=-a*(3*(1+Nb))/2;
% c=-3*a-2*b;
% d=-a-b-c;
% bl_=a*lam.^3+b*lam.^2+c*lam+d;
% bl(1)=0;
% bl(2:Nb+1)=bl_;
% bl(Nb+2)=1;
% s_C1=s_C1_;
% for i=1:Nb+2
%     s_C1(i)=r1(i)+bl(i)*(s_C1_(i)-r1(i));
%     s_C1(end+1-i)=r2(i)+bl(i)*(s_C1_(end+1-i)-r2(i));
% end
% 
% if plotting=='t'
%    plot(deriv(s_C1',1),'k') 
% %        plot(s_C1','k')
% %        plot(r1','b')
% 
% end

%% just checking if a b c d are correct 
% clearvars a b c d eqn1 eqn2 eqn3 eqn4
% syms a b c d
% eqn1 = a*lam(1).^3+b*lam(1).^2+c*lam(1)+d == 0
% eqn2 = a*lam(end).^3+b*lam(end).^2+c*lam(end)+d == 1
% eqn3 = 3*a*lam(1).^2+2*b*lam(1).^1+c == 0
% eqn4 = 3*a*lam(end).^2+2*b*lam(end).^1+c == 0
% [Ab,Bb] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4], [a,b,c,d]);
% Sol2 = mldivide(Ab,Bb);
% s_C2=double(Sol2(1,1))*lam.^3+double(Sol2(2,1))*lam.^2+double(Sol2(3,1))*lam+double(Sol2(4,1));
% Sol2
% figure
% plot(bl)
% hold on
% plot(s_C2)



end

