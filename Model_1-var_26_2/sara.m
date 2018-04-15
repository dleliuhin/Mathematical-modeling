function sara
clc 
clear all;
%Объявление переменных для блока А
syms p1 p2 p3 p4 p5 p6 x1 x2 ;
p1=0.5;
p3=1.0e+55;
p5=0.8;
p6=0;

dx1dt=-p1*x1+p2*(1-x1)*exp(x2/(x2/p3 + 1))
dx2dt=-p1*x2+p2*p4*(1-x1)*exp(x2/(x2/p3 + 1))-p5*(x2-p6)
x1m=x1;
x2m=x2;
p2m=p2;
p4m=p4;
%Блок аналитических преобразований
fprintf('Строим матрицу Якоби:');
jac=jacobian([dx1dt,dx2dt],[x1m,x2m])
fprintf('Строим сумму главных диагоналей:');
gld=simplify(jac(1,1)+jac(2,2))
SYSTEM=[dx1dt;dx2dt;gld]
fprintf('Умножим (1) уравнение на p4:');
SYSTEM(1)=SYSTEM(1)*p4
fprintf('Вычтем (1) уравнение из (2):');
act2=SYSTEM(2)-SYSTEM(1)
fprintf('Упростим получившееся уравнение');
act2=simplify(act2)
fprintf('Выразим x1:');
x1=solve(act2,x1)
fprintf('Подставим x1 в (1):');
act3=dx1dt;
act3=subs(act3,x1)
fprintf('Решим уравнение по параметру p2:');
p2=solve(act3,p2m)
fprintf('Упростим получившееся уравнение');
p2=simplify(p2)
fprintf('Подставим p2 и x1 в (3):');
act4=simplify(subs(SYSTEM(3),[p2m,x1m],[p2,x1]))
fprintf('Выразим p4 из (3):');
p4=solve(act4,p4)

P2=[];P4=[];L1=[];L2=[];J=[];

for x2t=1.0:0.1:7.0
    p4t=double(subs(p4,[x2m],[x2t]));
    p2t=double(subs(p2,[x2m,p4m],[x2t,p4t]));
    x1t=double(subs(x1,[x2m,p4m],[x2t,p4t]));
    jact=double(subs(jac,[x1m,x2m,p2m,p4m],[x1t,x2t,p2t,p4t]));
    dtjac=double(det(jact));
    lajac=double(eig(jact));    
    gldt=jact(1,1)+jact(2,2);
    if (p2t>0 && p4t>0)
        fprintf('\tx2=%1.1f | x1=%5.5f | p2=%5.5f | p4=%5.5f | |J|=%5.5f | G=%5.5f | L1=%s | L2=%s\n',...
                   x2t,x1t,p2t,p4t,dtjac,gldt,num2str(lajac(1)),num2str(lajac(2)));
    P4(end+1)=p4t;
    P2(end+1)=p2t;
    L1(end+1)=lajac(1);
    L2(end+1)=lajac(2);
    J(end+1)=dtjac;
    else
        fprintf('\tx2=%1.1f | x1=%5.5f | p2=%5.5f | p4=%5.5f | |J|=%5.5f | G=%5.5f | L1=%s | L2=%s\n',...
                   x2t,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
    end
    fprintf('\n');
end
TEMP2=[];TEMP4=[];count1=0;count2=0;
for k=1:length(P2)-1
    
    TEMP2=[P2(k),P2(k+1)];
    TEMP4=[P4(k),P4(k+1)];
    plot(TEMP2,TEMP4,'b.-');
    hold on;
    TEMP2=[];
    TEMP4=[];    
    if (imag(L1(k))<0 && imag(L1(k+1))>0) || (imag(L2(k))<0 && imag(L2(k+1))>0) || (J(k)<0 && J(k+1)>0)
    %С - на +
    TEMP2=[P2(k),P2(k+1)];
    TEMP4=[P4(k),P4(k+1)];
    plot(TEMP2,TEMP4,'r','LineWidth',4);
    TEMP2=[];
    TEMP4=[];
    hold on;    
    count1=count1+1; 
    end
    if (imag(L1(k))>0 && imag(L1(k+1))<0) || (imag(L2(k))>0 && imag(L2(k+1))<0) || (J(k)>0 && J(k+1)<0)
    %С + на -
    TEMP2=[P2(k),P2(k+1)];
    TEMP4=[P4(k),P4(k+1)];
    plot(TEMP2,TEMP4,'r','LineWidth',4);
    TEMP2=[];
    TEMP4=[];
    hold on;    
    count2=count2+1;    
    end
    title('(p2,p4)','fontsize',15);
    xlabel('p2','fontsize',15); 
    ylabel('p4 ','fontsize',15,'rotation',0); 
end
fprintf('Value of points Andronov-Hoppf:  %i',count1+count2);
P2=[];P4=[];L1=[];L2=[];
end        