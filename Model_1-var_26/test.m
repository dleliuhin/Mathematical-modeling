function test(bx1,bx2,bp2,bp4,bl1,bl2,bj)
%Объявление переменных
p1=0.5;
p3=1.0e+55;
p5=0.8;
p6=0;
syms x1 x2 p2 p4;
dx1dt=-p1*x1+p2*(1-x1)*exp(x2/(x2/p3 + 1));
dx2dt=-p1*x2+p2*p4*(1-x1)*exp(x2/(x2/p3 + 1))-p5*(x2-p6);
maj=jacobian([dx1dt,dx2dt],[x1,x2]);
gld=simplify(maj(1,1)+maj(2,2));
laj=eig(maj);
l1=laj(1);
l2=laj(2);
fprintf('\n\t БЛОК ПРОВЕРКИ\n')
%Блок проверки
for b=1:length(bp4)
    f1=double(subs(dx1dt,[x1,x2,p2],[bx1(b),bx2(b),bp2(b)]));
    f2=double(subs(dx2dt,[x1,x2,p2,p4],[bx1(b),bx2(b),bp2(b),bp4(b)]));
    J=double(subs(maj,[x1,x2,p2,p4],[bx1(b),bx2(b),bp2(b),bp4(b)]));
    DJ=det(J);   
    G=double(subs(gld,[x1,x2,p2,p4],[bx1(b),bx2(b),bp2(b),bp4(b)]));
fprintf('\t x2=%4.1f | x1=%4.1f | p2=%5.5f | p4=%5.5f | |J|=%5.5f | G=%5.5f | dx1dt=%5.5f | dx2dt=%5.5f \n', ...
                  bx2(b),bx1(b),bp2(b),bp4(b),DJ,G,f1,f2);        
end    
end