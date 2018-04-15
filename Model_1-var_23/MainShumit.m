function M=MainShumit
%������� ���� ������ � ������� ����������
clc
clear all;
%���������� ����������
syms p1 p3 p4 p5 p6;
p1=1;
p3=20;
p4=10;
p5=0.6;
p6=-5;

syms x1 x2 p2;
syms maj majt;
sym la;
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
fprintf('������ ������� ����� �� ��������� p2:');
jacp2=jacobian([dx1dt,dx2dt],[p2m])
fprintf('������ ������������ ������� �����:');
deter=det(jac)
fprintf('��������� x1 � (3):');
act4=subs(deter,x1)
fprintf('����� ��������� �� ��������� p2:');
nextp2=solve(act4,p2m)
fprintf('�������� ������������ ���������');
nextp2=simplify(nextp2)
fprintf('��������� x1 � p2 � (1):');
tmp=subs(SYSTEM2(1),[x1m,p2m],[x1,nextp2])
fprintf('����� ��������� �� ���������� x2:');
tmp1x2=double(solve(tmp,x2m))
fprintf('��������� x1 � p2 � (3)-������������ ������� ����� (deter):');
tmp=subs(deter,[x1m,p2m],[x1,nextp2])
fprintf('����� ��������� �� ���������� x2:');
tmp2x2=double(solve(tmp,x2m))
fprintf('������� ��� ������������ �������� �2 � ��������� ������:');
X2=[];X2=[tmp1x2',tmp2x2]
fprintf('����������� �� �������� ������� �2:');
X2=sort(X2)

P2=[];TX2=[];X1=[];LA1=[];LA2=[];
JP2=[jacp2(1) jac(1,2);jacp2(2) jac(2,2)];

fprintf('\n\t������� �����:\n\n');
for k=1:length(X2)
   x1t=double(subs(x1,[x2m],[X2(k)]));
   try
   p2t=double(subs(p2,[x2m],[X2(k)]));
   catch
   end    
   jact=double(subs(jac,[x1m,x2m,p2m],[x1t,X2(k),p2t]));
   detj=double(det(jact));
   jacp=double(subs(JP2,[x1m,x2m,p2m],[x1t,X2(k),p2t]));
   detjp=double(det(jacp));   
   gld=double(jact(1,1)+jact(2,2));
   la=double(eig(jact));
   f1=double(subs(dx1dt,[x1m,x2m,p2m],[x1t,X2(k),p2t]));
   f2=double(subs(dx2dt,[x1m,x2m,p2m],[x1t,X2(k),p2t]));   
   if (p2t>0 && x1t>0 && x1t<1)
       fprintf('\tx2=%5.5f | x1=%5.5f | p2=%5.5f | L1=%s | L2=%s | |J|=%5.5f | dx1dt=%5.5f | dx2dt=%5.5f\n', ...
                  X2(k),x1t,p2t,num2str(la(1)),num2str(la(2)),detj,f1,f2);
       TX2(end+1)=X2(k);
       X1(end+1)=x1t;
       P2(end+1)=p2t;
       LA1(end+1)=la(1);
       LA2(end+1)=la(2);
   else
       fprintf('\tx2=%5.5f | x1=%5.5f | p2=%5.5f | L1=%s | L2=%s | |J|=%5.5f | G=%5.5f \n', ...
                  X2(k),NaN,NaN,NaN,NaN,NaN,NaN);
   end
   
   fprintf('\n');
end

%���� ��������
fprintf('\t���� ��������: \n');
for q=1:length(X2)
    %����������� � ������ ��������� ��� ��������� ��������
    f1=double(subs(dx1dt,[p2m,x1m,x2m],[P2(q),X1(q),TX2(q)]));
    %����������� �� ������ ��������� ��� ��������� ��������
    f2=double(subs(dx2dt,[p2m,x1m,x2m],[P2(q),X1(q),TX2(q)]));
    %����������� � ������� ����� ��� ��������� ��������
    DJ=double(subs(jac,[p2m,x1m,x2m],[P2(q),X1(q),TX2(q)]));
    %��������� ������������ ������� �����
    dt=det(DJ);
    %������� ���� �������� � ��������� ����
fprintf('\t x2=%4.1f | x1=%5.5f | p2=%5.5f | |J|=%5.5f | dx1dt=%5.5f | dx2dt=%5.5f \n', ...
                  TX2(q),X1(q),P2(q),dt,f1,f2);    
end
end