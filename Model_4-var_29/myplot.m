function [fig]=myplot(p4t,fig,X,P1,L1,L2)
kolv=0;%������� ����������� � ������� ������������ ����������
kolk=0;%������� ����������� � ������� ���������� ���������-�����
stac=0;%������� ������������ ���������� �����
nstac=0;%������� ������������ ������������ �����
ost=0;%������� �����, �� ������������� �� ������ �������

%� ����������� �� �������� ������� ������ ��������� ���������� strx
%�������� 1 ��� 2
if (rem(fig,2)==1)
    strx=1;
elseif (rem(fig,2)==0)
    strx=2;
end

%������ ������ � ����� ���� � �������� fig
figure(fig);

%������� �������� ������� ��� ���������� �����������
%����� ����� ������� �� �������
TEMPX=[];TEMP1=[];

%���� � �������� ���������� ������ �� ������� ������� length(�1)-1
%�� 1 �� ����� ������ �������, �������� �1
%(����� ���� �������� ����������)
for k=1:length(P1)-1
    
    %������ ��������� ��� ����������� ����� �� �������
    TEMPX=[X(k),X(k+1)];
    TEMP1=[P1(k),P1(k+1)];
    %������������ ������� ������ ������� �����
    plot(TEMP1,TEMPX,'k-');
    %������� hold on �������� ����� ���������� �������� ������� 
    hold on;
    %����������� �������� ������� �� ��������
    TEMP2=[];TEMP4=[];    
    
%.....������ ������� ��� 4 ����� ����� � ���������� �� �� ��������........
    
    if (real(L1(k))<0 && real(L2(k))<0)
        %������������ ���������� �����
        %������ "���������" �������� �����
        plot(P1(k),X(k),'Marker','p','MarkerEdgeColor','g', ...
                                     'MarkerFaceColor','g');
        %������� hold on �������� ����� ���������� �������� �������                                  
        hold on;
        stac=stac+1;%����������� ������� �����
        
    elseif ((real(L1(k))>0 && real(L2(k))>0) || ...
            (real(L1(k))<0 && real(L2(k))>0) || ...
            (real(L1(k))>0 && real(L2(k))<0))
        %������������ ������������ �����
        %������ "o" �������� �����    
        plot(P1(k),X(k),'Marker','o','MarkerEdgeColor','r', ...
                                            'MarkerFaceColor','r');                                        
        %������� hold on �������� ����� ���������� �������� �������                                          
        hold on;    
        nstac=nstac+1;%����������� ������� �����            
    else
        %������ �����
        %������ "�" ������� �����
        plot((P1(k)),(X(k))','ko');
        %������� hold on �������� ����� ���������� �������� �������  
        hold on;
        ost=ost+1;%����������� ������� �����           
    end
    
    
    	if (real(L1(k))<0 && real(L1(k+1))>0) || ...
           (real(L2(k))<0 && real(L2(k+1))>0) 
        %���������� � ������� ������������ ���������� (�������� � - �� +)
        %������ "������" ������ �����
        
       %���������� � �������� ������� �������� ����������� ����� �� �������
            TEMPX=[X(k),X(k+1)];
            TEMP1=[P1(k),P1(k+1)];
            %������ ���������� ����� ������ ����� ������� 2 ��.
            plot(TEMP1,TEMPX,'Color','b','LineWidth',2.0,'LineStyle','--');
            %������ ������ ���������� ����� ��� �����������
            plot((P1(k)+P1(k+1))/2,(X(k)+X(k+1))/2, ...
                'Marker','h','MarkerSize',11,'MarkerEdgeColor','b', ...
                                             'MarkerFaceColor','b'); 
            %������� hold on �������� ����� ���������� �������� �������  
            hold on;
            kolv=kolv+1;%����������� ������� �����
            %����������� �������� ������� �� ��������
            TEMPX=[];TEMP1=[];
        end
        
        if (real(L1(k))>0 && real(L1(k+1))<0) || ...
           (real(L2(k))>0 && real(L2(k+1))<0)
        %���������� � ������� ������������ ���������� (�������� � + �� -)
        %������ "������" ������ �����            
        
        %���������� � �������� ������� �������� ����������� ����� �� �������
            TEMPX=[X(k),X(k+1)];
            TEMP1=[P1(k),P1(k+1)];
            %������ ���������� ����� ������ ����� ������� 2 ��.
            plot(TEMP1,TEMPX,'Color','b','LineWidth',2.0,'LineStyle','--');
            %������ ������ ���������� ����� ��� �����������
            plot((P1(k)+P1(k+1))/2,(X(k)+X(k+1))/2, ...
                'Marker','h','MarkerSize',11,'MarkerEdgeColor','b', ...
                                             'MarkerFaceColor','b');     
            %������� hold on �������� ����� ���������� �������� �������  
            hold on;    
            kolv=kolv+1;%����������� ������� �����
            %����������� �������� ������� �� ��������
            TEMPX=[];TEMP1=[];
        end

    	if (imag(L1(k))<0 && imag(L1(k+1))>0) || ...
           (imag(L2(k))<0 && imag(L2(k+1))>0) 
        %���������� � ������� ����������� ���������� ���������-����� (�������� � - �� +)
        %������ "�������" ����� ����                
        
       %���������� � �������� ������� �������� ����������� ����� �� �������
            TEMPX=[X(k),X(k+1)];
            TEMP1=[P1(k),P1(k+1)];
            %������ ���������� ����� ����� ���� ������� 2 ��.
            plot(TEMP1,TEMPX,'Color','�','LineWidth',2.0,'LineStyle','--');
            %������ ������ ���������� ����� ��� �����������
            plot((P1(k)+P1(k+1))/2,(X(k)+X(k+1))/2, ...
                'Marker','s','MarkerSize',11,'MarkerEdgeColor','�', ...
                                             'MarkerFaceColor','�');     
            %������� hold on �������� ����� ���������� �������� �������  
            hold on;
            kolk=kolk+1;%����������� ������� �����
            %����������� �������� ������� �� ��������
            TEMPX=[];TEMP1=[];
        end
        if (imag(L1(k))>0 && imag(L1(k+1))<0) || ...
           (imag(L2(k))>0 && imag(L2(k+1))<0)
        %���������� � ������� ����������� ���������� ���������-����� (�������� � + �� -)
        %������ "�������" ����� ����         
        
       %���������� � �������� ������� �������� ����������� ����� �� �������
            TEMPX=[X(k),X(k+1)];
            TEMP1=[P1(k),P1(k+1)];
            %������ ���������� ����� ����� ���� ������� 2 ��.
            plot(TEMP1,TEMPX,'Color','�','LineWidth',2.0,'LineStyle','--');
            %������ ������ ���������� ����� ��� �����������
            plot((P1(k)+P1(k+1))/2,(X(k)+X(k+1))/2, ...
                'Marker','s','MarkerSize',11,'MarkerEdgeColor','�', ...
                                             'MarkerFaceColor','�');     
            %������� hold on �������� ����� ���������� �������� �������  
            hold on;    
            kolk=kolk+1;%����������� ������� �����
            %����������� �������� ������� �� ��������
            TEMPX=[];TEMP1=[];
        end        
end
%������� grid on ������� ������������ ����� �� ������� ���    
grid on;
%��������� ����� ��� �������� ��� ������ title
title(sprintf('������ ����������� x%i (p1) ��� p4=%1.1f:',strx,p4t),'fontsize',15);
%�������� ����� �� ��� x � tau
xlabel('p1','fontsize',15); ylabel(sprintf('x%i ',strx),'fontsize',15,'rotation',0);
%������� ����������� ����� ������ �����    
fprintf('\n\t��� p4=%1.1f',p4t);
fprintf('\n\t����������� ������������ ���������� ����� ��� x%i(p1):%i',strx,stac);
fprintf('\n\t����������� ������������ ������������ ����� ��� x%i(p1):%i',strx,nstac);
fprintf('\n\t����������� ����� ������������ ���������� ��� x%i(p1):%i',strx,kolv);
fprintf('\n\t����������� ����� ���������� ���������-����� ��� x%i(p1):%i',strx,kolk);
fprintf('\n\t����������� ������ ����� ��� x%i(p1):%i',strx,ost);
fprintf('\n');
fig=fig+1;%����������� ������� fig
end