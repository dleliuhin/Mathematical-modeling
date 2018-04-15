function Test(pra,prtau,prx1,prx2)
%Объявление переменных
p1=1;
p3=20;
p4=10;
p6=-5;
k0=1;
syms a x1 x2 tau;
syms maj majt;
sym ls;
p2=k0*tau;
p5=a*tau;
dx1dt=-p1*x1+p2*(1-x1)*exp(x2/(1+x2/p3)); 
dx2dt=-p1*x2+p2*p4*(1-x1)*exp(x2/(1+x2/p3))-p5*(x2-p6);
maj=jacobian([dx1dt,dx2dt],[x1,x2]);
fprintf('\n\t БЛОК ПРОВЕРКИ\n')
%Блок проверки
for b=1:length(prtau)
    f1=double(subs(dx1dt,[tau,x1,x2],[prtau(b),prx1(b),prx2(b)]));
    f2=double(subs(dx2dt,[a,tau,x1,x2],[pra(b),prtau(b),prx1(b),prx2(b)]));
    J=double(subs(maj,[a,x1,x2,tau],[pra(b),prx1(b),prx2(b),prtau(b)]));
    DJ=det(J);
fprintf('\ta=%i x2=%4.1f | x1=%5.5f | tau=%5.5f | |J|=%5.5f | dx1dt=%5.5f | dx2dt=%5.5f \n', ...
                  pra(b),prx2(b),prx1(b),prtau(b),DJ,f1,f2);    
end
end