function myplot(P2,P4,L1,L2,J)
%�������������� �������� ������� ��� ��������� ���� ����� ������� ��������
%�������� ��� �2 � ��� �4
TEMP2=[];TEMP4=[];
%������� ����������-�������� ����� ���������� ���������-�����
%����� ����� ���������� �.-�. ��� ��������� ����� ������ ����� �����������
%�������� ������� ����� ��� ������������ ������� ����� � - �� +  
count1=0;
%����� ����� ���������� �.-�. ��� ��������� ����� ������ ����� �����������
%�������� ������� ����� ��� ������������ ������� ����� � + �� -  
count2=0;
%�������� ���� �� 1 �� ����� ������� �������� ���������� -1
%..........................................................................
for k=1:length(P2)-1
    
    %������� ��� �����, ���������� �� �������
    %���������� � �������� ������ ��� ����� ������� �������� �2
    TEMP2=[P2(k),P2(k+1)];
    %���������� � �������� ������ ��� ����� ������� �������� �4
    TEMP4=[P4(k),P4(k+1)];
    %������ ������ �������(�� ���� �� ���� �������� �2 � �4), �.�. �� 2 ����� �� ���
    %��� ����������� ������� ������ ���������� ����� ������
    plot(TEMP2,TEMP4,'b.-');
    %������� hold on - ������ ����� �������� ������ ������
    hold on;
    %�������� �������� ������� �������� �2 � �4
    TEMP2=[];TEMP4=[];    
    
%..........................................................................    
    %��������� ��� ����� ������� ����� �� ���������� ������� - �� +
    if (imag(L1(k))<0 && imag(L1(k+1))>0) || (imag(L2(k))<0 && imag(L2(k+1))>0) || (J(k)<0 && J(k+1)>0)
        %���������� � �������� ������ ��� ����� ������� �������� �2
        TEMP2=[P2(k),P2(k+1)];
        %���������� � �������� ������ ��� ����� ������� �������� �4
        TEMP4=[P4(k),P4(k+1)];
        %������ ������ �������(�� ���� �� ���� �������� �2 � �4), �.�. �� 2 ����� �� ���
        %��� ����������� ������� ������ ���������� ������ ������� ������
        plot(TEMP2,TEMP4,'r','LineWidth',5);
        %������� hold on - ������ ����� �������� ������ ������
        hold on;    
        %�������� �������� ������� �������� �2 � �4
        TEMP2=[];TEMP4=[];
        %���� ������� �����������, ����������� ������� �� 1
        count1=count1+1; 
    end
    
%..........................................................................    
    %��������� ��� ����� ������� ����� �� ���������� ������� + �� -
    if (imag(L1(k))>0 && imag(L1(k+1))<0) || (imag(L2(k))>0 && imag(L2(k+1))<0) || (J(k)>0 && J(k+1)<0)
        %���������� � �������� ������ ��� ����� ������� �������� �2
        TEMP2=[P2(k),P2(k+1)];
        %���������� � �������� ������ ��� ����� ������� �������� �4
        TEMP4=[P4(k),P4(k+1)];
        %������ ������ �������(�� ���� �� ���� �������� �2 � �4), �.�. �� 2 ����� �� ���
        %��� ����������� ������� ������ ���������� ������ ������� ������
        plot(TEMP2,TEMP4,'g','LineWidth',5);
        %������� hold on - ������ ����� �������� ������ ������
        hold on;    
        %�������� �������� ������� �������� �2 � �4
        TEMP2=[];TEMP4=[];
        %���� ������� �����������, ����������� ������� �� 1
        count2=count2+1;    
%..........................................................................        
    end
    %��������� ����� ��� �������� ��� ������ title
    title('(p2,p4)','fontsize',15);
    %�������� ��� x
    xlabel('p2','fontsize',15);
    %�������� ��� y
    ylabel('p4 ','fontsize',15,'rotation',0); 
end
%������� ���������� ����� ���������� ���������-�����
fprintf('���������� ����� ���������� ���������-�����:  %i \n',count1+count2);
end