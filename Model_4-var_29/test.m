function test(X1,X2,P1,P4)
%Объявление всех используемых переменных
syms p1 p2 p3 p4 p5;
syms x1 x2;
%Присваиваем некоторым переменным значения из условия задания
p2=1.5;
p3=3;
p5=0.01;
%Запиcываем два главных уравнения из задания
dx1dt=p1*(p5+x1^p3)/(1+x1^p3)-x1*(1+x2);
dx2dt=x1*(p2+x2)-p4*x2;
maj=jacobian([dx1dt,dx2dt],[x1,x2]);
gld=simplify(maj(1,1)+maj(2,2));

fprintf('\n\t БЛОК ПРОВЕРКИ\n')
for b=1:length(P1)
    f1=double(subs(dx1dt,[x1,x2,p1],[X1(b),X2(b),P1(b)]));
    f2=double(subs(dx2dt,[x1,x2,p4],[X1(b),X2(b),P4(b)]));
    J=double(subs(maj,[x1,x2,p1,p4],[X1(b),X2(b),P1(b),P4(b)]));
    DJ=det(J);   
    G=double(subs(gld,[x1,x2,p1,p4],[X1(b),X2(b),P1(b),P4(b)]));
fprintf('p4=%3.1f | x1=%3.2f | x2=%5.3f | p1=%5.3f | |J|=%5.3f | G=%5.3f | dx1dt=%5.3f | dx2dt=%5.3f \n', ...
                  P4(b),X1(b),X2(b),P1(b),DJ,G,f1,f2);        
end
end