function M=Main
%������� ���� ������ � ������� ����������
clc
clear all;
%���������� ����������
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

x1j=x1;%����������� �������� x1 � �������� ����������
x2j=x2;%����������� �������� x2 � �������� ����������
%������� ���� OUTValues.dat ��� ���������� ������ ������� ������ � ����
outfile=fopen('D:\���.������\������\OUTValues.dat','w');
%--------------------------------------------------------------------------
%���� ������������� ��������������
fprintf('������� (1) ��������� �� p4:');
temp1=dx1dt*p4
fprintf('������ (1) ��������� �� (2):');
temp2=dx2dt-temp1
fprintf('������� x1:');
x1=solve(temp2,x1)
fprintf('��������� x1 � (1):');
temp3=dx1dt;
temp3=subs(temp3,x1)
fprintf('����� ���������� ��������� �� ���������� tau:');
taut=solve(temp3,tau)
fprintf('������ ������� �����:');
maj=jacobian([dx1dt,dx2dt],[x1j,x2j])
%--------------------------------------------------------------------------
%�������� ����� (������,����, ����� ����������� � ��������)
%������� �������� ������ 1*2 ��� ������ ������ ��������� (tau) � ����
kor=[0 0];
%������� �������� ������ 2*2 ��� ������ ������� ����� � ����
majt=[0 0;0 0];
%������� ������ ������� ��� ������ ���� ������ � ���� �������
KP=[];A=[];X2=[];TAU=[];K=[];X1=[];LAM1=[];LAM2=[];
%������� ������ �������� ������� ��� ������ �������� �� ������
BUFA=[];BUFTAU=[];BUFX1=[];BUFX2=[];BUFLAM1=[];BUFLAM2=[];
GRAFA=[];GRAFTAU=[];GRAFX1=[];GRAFX2=[];GRAFLAM1=[];GRAFLAM2=[];
%������� ������ �������� ������� ��� ���������� ����� ��������
BUFTAU1=[];BUFTAU2=[];BUFX11=[];BUFX12=[];BUFX21=[];BUFX22=[];
TEMPLAM11=[];TEMPLAM12=[];TEMPLAM21=[];TEMPLAM22=[];BUFA1=[];BUFA2=[];
iter=1;%������� ������ ����������� �����
fig=1;%����������, ���������� �� ����� ������ ����������� �������
at=1;%���������� ���������� �, �������� �������� � ����� �=1 2 4
%������ �����, ���� � �� ������ ������ 4
while at<=4
    for x2t=-5.0:0.1:10.0
        try
        kor=double(subs(taut,[x2,a],[x2t,at]));
        catch
        %���� try-catch ��������� ��� ���������� ������ � ��������������
        %���� '������� �� ����' ��� ������� ����������� ��������� �� tau
        end
        %����� ��������� �������� �� ������ tau ����������� ������
        if imag(kor)~=0
            kor=abs(kor);
        end
        %���� ����������� ���� ����� �������� ���������� � ���� ������ 
        fprintf('a=%d x2=%4.1f\n',at,x2t);
        %��� ��� � ��� �� 2 �������� tau=>���� �� kor(k) (��� tau)
        for k=1:length(kor) %����� kor ��������� �������� �� 1 �� 2
            %����������� �������� a,x2,tau � �������� ���������� x1t
            x1t=double(subs(x1,[a,x2,tau],[at,x2t,kor(k)]));
            %����������� �������� a,x1,x2,tau � �������� ������� ����� majt
            majt=double(subs(maj,[a,x1j,x2j,tau],[at,x1t,x2t,kor(k)]));
            %��������� ����� ��������� ������� ��������� ������� majt
            sd=majt(1,1)+majt(2,2);
            %������� ������ ����������� �������� ������� ����� majt
            ls=double(eig(majt));
            %��� ������ ������ � ����� ������� ��������� ���������� ��� 
            %�������� � ������� 
            A(end+1)=at;
            X2(end+1)=x2t;
            KP(end+1)=iter;
            %����������� ������� ������ ����������� �����
            iter=iter+1;
            %������ ������� �� �������
            if (kor(k)>0 && x1t>0 && x1t<1)
            %���� ����������� ���� ����� �������� ���������� � ���� ������ 
                 fprintf('\ttau(%d)= %9.5f | x1=%9.5f | L1=%s | L2=%s \n', ...
                 k,kor(k),x1t,num2str(ls(1)),num2str(ls(2)));
            %����� ��������� ��������� ����������� �������� � ��������
             if (k==1)
                TAU(end+1)=kor(k);                
                K(end+1)=k;
                X1(end+1)=x1t;
                LAM1(end+1)=ls(1);
                LAM2(end+1)=ls(2);
                
                BUFA1(end+1)=at;                
                BUFTAU1(end+1)=kor(k);
                BUFX11(end+1)=x1t;
                BUFX21(end+1)=x2t;
                TEMPLAM11(end+1)=ls(1);
                TEMPLAM12(end+1)=ls(2);
            elseif (k==2)
                TAU(end+1)=kor(k);
                K(end+1)=k;
                X1(end+1)=x1t;
                LAM1(end+1)=ls(1);
                LAM2(end+1)=ls(2);

                BUFA2(end+1)=at;                
                BUFTAU2(end+1)=kor(k);
                BUFX12(end+1)=x1t;
                BUFX22(end+1)=x2t;
                TEMPLAM21(end+1)=ls(1);
                TEMPLAM22(end+1)=ls(2);
             end
            %���� �� ����� tau �����������, ����������� ������ ��������
            else
                 TAU(end+1)=NaN;
                 K(end+1)=NaN;
                 X1(end+1)=NaN;
                 LAM1(end+1)=NaN;
                 LAM2(end+1)=NaN;                 
                 fprintf('\t%d-� ������ tau �����������\n',k);
            continue
            %���� k-�� ������ �����������, �������� ������� ���-�� �����
                 iter=iter-1;
            end
        end
        fprintf('\n');
    end 
    %���� �������� "���-��������������" ����� �� ��������:
    %�������������� ����� ���������� ������ ������ ������ ���������
    %� ������������� ��� �������� ������� ��������
    BUFA2=fliplr(BUFA2);
    BUFTAU2=fliplr(BUFTAU2);
    BUFX12=fliplr(BUFX12);
    BUFX22=fliplr(BUFX22);
    TEMPLAM21=fliplr(TEMPLAM21);
    TEMPLAM22=fliplr(TEMPLAM22);
    %� ��������� � ��������� ������ ������ ���������
    BUFA=[BUFA1,BUFA2];
    BUFTAU=[BUFTAU1,BUFTAU2];
    BUFX1=[BUFX11,BUFX12];
    BUFX2=[BUFX21,BUFX22];
    BUFLAM1=[TEMPLAM11,TEMPLAM21];
    BUFLAM2=[TEMPLAM12,TEMPLAM22];
    %�������� �������� ������ �������� ��� ������ ����� ���� ��������
    GRAFA=[GRAFA,BUFA];
    GRAFTAU=[GRAFTAU,BUFTAU];
    GRAFX1=[GRAFX1,BUFX1];
    GRAFX2=[GRAFX2,BUFX2];
    GRAFLAM1=[GRAFLAM1,BUFLAM1];
    GRAFLAM2=[GRAFLAM2,BUFLAM2];
    %���������� ����������� ������� ���������� �������� ��� x1 � x2
    fig=mylogplot(at,fig,BUFTAU,BUFLAM1,BUFLAM2,BUFX1);
    fig=mylogplot(at,fig,BUFTAU,BUFLAM1,BUFLAM2,BUFX2);
    %�������� �������, ����� ����� ��� ������ �������� �
    %�� ������ � ����� ���������
    BUFA=[];BUFTAU=[];BUFX1=[];BUFX2=[];BUFLAM1=[];BUFLAM2=[]; 
    BUFTAU1=[];BUFTAU2=[];BUFX11=[];BUFX12=[];BUFX21=[];BUFX22=[];
    TEMPLAM11=[];TEMPLAM12=[];TEMPLAM21=[];TEMPLAM22=[];BUFA1=[];BUFA2=[];
    %�������� ����������� � �� 2
     at=at*2;
end
%�������������� ������� �������� ���� ����� � ������� � ���� ������ 
T=table;
T.N=KP';T.a=A';T.x2=X2';T.tau=TAU';
T.k=K';T.x1=X1';T.LAM1=LAM1';T.LAM2=LAM2'
%�������������� ������� �������� ����� �� ������� ������ �������
G=table;
G.N=(1:length(GRAFTAU))';G.a=GRAFA';G.x2=GRAFX2';G.tau=GRAFTAU';
G.x1=GRAFX1';G.LAM1=(double(GRAFLAM1))';G.LAM2=(single(GRAFLAM2))'
%������ � ��������� ���� ������� �����, ����� ��� ���� ��������: 
figure(fig);
imshow('D:\���.������\������\Legend.png');
%���������� ������� � ��������� ���� OUTValues.dat 
writetable(T,'D:\���.������\������\OUTValues.dat');
%��������� ����
fclose(outfile);
%���� ��������
Test(GRAFA,GRAFTAU,GRAFX1,GRAFX2);
end