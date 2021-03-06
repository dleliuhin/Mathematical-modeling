function M=AMainStarih
%������� ���� ������ � ������� ����������
clc
clear all;
%���������� ���������� ��� ����� �
syms p1 p3 p4 p5 p6;
p1=1;
p3=20;
p4=10;
p5=0.6;
p6=-5;
%����� ��������� ����������� ���������� 
syms x1 x2 p2;
syms maj majt;
sym la;
%�������� ��������� 
dx1dt=-p1*x1+p2*(1-x1)*exp(x2/(1+x2/p3))
dx2dt=-p1*x2+p2*p4*(1-x1)*exp(x2/(1+x2/p3))-p5*(x2-p6)
x1m=x1;
x2m=x2;
p2m=p2;
%���� ������������� ��������������
SYSTEM2=[dx1dt;dx2dt]
fprintf('������� (1) ��������� �� p4:');
SYSTEM2(1)=SYSTEM2(1)*p4
fprintf('������ (1) ��������� �� (2):');
act2=SYSTEM2(2)-SYSTEM2(1)
fprintf('�������� ������������ ���������');
act2=simplify(act2)
fprintf('������� x1:');
x1=solve(act2,x1)
fprintf('��������� x1 � (1):');
act3=dx1dt;
act3=subs(act3,x1)
fprintf('����� ��������� �� ��������� p2:');
p2=solve(act3,p2)
fprintf('�������� ������������ ���������');
p2=simplify(p2)
fprintf('������ ������� �����:');
jac=jacobian([dx1dt,dx2dt],[x1m,x2m])
%�������������� �������, � ������� ����� ���������� ���������� ����������
P2=[];X2=[];X1=[];LA1=[];LA2=[];
%�������� ����������� �2 �� ��������� ���������
for x2t=-2.0:0.1:5.0
   %����������� �2 � ������� ���������� �1 
   x1t=double(subs(x1,[x2],[x2t]));
   try
   %����������� �2 � ������� ��������� �2
   p2t=double(subs(p2,[x2],[x2t]));
   catch
   end      
   %����������� ����������� �������� �1,�2,�2 � ������� �����
   jact=double(subs(jac,[x1m,x2m,p2m],[x1t,x2t,p2t]));
   %��������� ����������� �������� ������� �����
   la=double(eig(jact));
   %������ ���������� ������� ������ ��� ���������� "������" �����
   if (p2t>0 && x1t>0 && x1t<1)
       %������� ��� ���������� ����� � ��������� ����
       fprintf('\tx2=%4.1f | x1=%5.5f | p2=%5.5f | L1=%s | L2=%s \n', ...
                  x2t,x1t,p2t,num2str(la(1)),num2str(la(2)));
       %��������� �������� � �������
       X2(end+1)=x2t;
       X1(end+1)=x1t;
       P2(end+1)=p2t;
       LA1(end+1)=la(1);
       LA2(end+1)=la(2);
   else
       %���� ������� �� ����������� ������� NaN
       fprintf('\tx2=%4.1f | x1=%5.5f | p2=%5.5f | L1=%s | L2=%s \n', ...
                  x2t,NaN,NaN,NaN,NaN);
   end
   %fprintf('\n');
end
%���������� str ������ ��� ���������� ����������� �������� � ������ ������
str=1;
%������ ������ �1(�2)
myplotA(P2,X1,LA1,LA2,str);
str=2;
%������ ������ �2(�2)
myplotA(P2,X2,LA1,LA2,str);
%���� ��������
fprintf('\t���� ��������: \n');
for q=1:length(P2)
    %����������� � ������ ��������� ��� ��������� ��������
    f1=double(subs(dx1dt,[p2m,x1m,x2m],[P2(q),X1(q),X2(q)]));
    %����������� �� ������ ��������� ��� ��������� ��������
    f2=double(subs(dx2dt,[p2m,x1m,x2m],[P2(q),X1(q),X2(q)]));
    %����������� � ������� ����� ��� ��������� ��������
    DETMAJ=double(subs(jac,[p2m,x1m,x2m],[P2(q),X1(q),X2(q)]));
    %��������� ������������ ������� �����
    deter=det(DETMAJ);
    %������� ���� �������� � ��������� ����
fprintf('\t x2=%4.1f | x1=%5.5f | p2=%5.5f | |J|=%5.5f | dx1dt=%5.5f | dx2dt=%5.5f \n', ...
                  X2(q),X1(q),P2(q),deter,f1,f2);
end
end