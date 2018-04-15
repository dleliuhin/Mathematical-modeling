function [fig]=mylogplot(at,fig,BUFTAU,BUFLAM1,BUFLAM2,BUFX)
kolv=0;%������� ����� ������������ ����������
kolk=0;%������� ����� ���������� ���������-�����
stac=0;%������� ������������ ���������� �����
nstac=0;%������� ������������ ������������ �����
ost=0;%������� �����, �� �������� �� � ���� �������
%� ����������� �� �������� ������� ������ ��������� ���������� strx
%�������� 1 ��� 2
if (rem(fig,2)==1)
    strx=1;
elseif (rem(fig,2)==0)
    strx=2;
end
    %������ ������ � ����� ���� � �������� fig
    figure(fig);
    %���� � �������� ���������� ������ �� ������� ������� length(BUFTAU)-1
    for g=1:length(BUFTAU)-1
    %������ ������� ��� 4 ����� ����� � ���������� �� �� ��������
    %� ��������������� �������� �� ��� tau
        if (real(BUFLAM1(g))<0 && real(BUFLAM2(g))<0)
            %������������ ���������� �����
            %������ "���������" �������� �����
            semilogx(BUFTAU(g),BUFX(g),'Marker','p','MarkerEdgeColor','g', ...
                                                    'MarkerFaceColor','g');
            %������� hold on �������� ����� ���������� �������� ������� 
            hold on;
            stac=stac+1;%����������� ������� �����
%--------------------------------------------------------------------------            
        elseif ((real(BUFLAM1(g))>0 && real(BUFLAM2(g))>0) || ...
                (real(BUFLAM1(g))<0 && real(BUFLAM2(g))>0) || ...
                (real(BUFLAM1(g))>0 && real(BUFLAM2(g))<0))
            %������������ ������������ �����
            %������ "o" �������� �����            
            semilogx(BUFTAU(g),BUFX(g),'Marker','o','MarkerEdgeColor','r', ...
                                                    'MarkerFaceColor','r');
            hold on;    
            nstac=nstac+1;%����������� ������� �����
%--------------------------------------------------------------------------            
        else
           %������ �����
           %������ "�" ������� �����
           semilogx((BUFTAU(g))',(BUFX(g))','ko');
           hold on;
           ost=ost+1;%����������� ������� �����
        end
%----------------------------------------------------------------------------------------------
		if (real(BUFLAM1(g))<0 && real(BUFLAM1(g+1))>0) || ...
                (real(BUFLAM2(g))<0 && real(BUFLAM2(g+1))>0) 
            %���������� � ������� ������������ ���������� (�������� � - �� +)
            %������ "������" ������ �����
            semilogx((BUFTAU(g)+BUFTAU(g+1))/2,(BUFX(g)+BUFX(g+1))/2, ...
                'Marker','h','MarkerSize',11,'MarkerEdgeColor','b', ...
                                             'MarkerFaceColor','b');
            hold on;
            kolv=kolv+1;%����������� ������� �����
        end
        if (real(BUFLAM1(g))>0 && real(BUFLAM1(g+1))<0) || ...
                (real(BUFLAM2(g))>0 && real(BUFLAM2(g+1))<0)
            %���������� � ������� ������������ ���������� (�������� � + �� -)
            %������ "������" ������ �����            
            semilogx((BUFTAU(g)+BUFTAU(g+1))/2,(BUFX(g)+BUFX(g+1))/2, ...
                'Marker','h','MarkerSize',11,'MarkerEdgeColor','b', ...
                                             'MarkerFaceColor','b');
            hold on;    
            kolv=kolv+1;%����������� ������� �����
        end
%--------------------------------------------------------------------------        
		if (imag(BUFLAM1(g))<0 && imag(BUFLAM1(g+1))>0) || ...
                (imag(BUFLAM2(g))<0 && imag(BUFLAM2(g+1))>0) 
            %���������� � ������� ����������� ���������� ���������-����� (�������� � - �� +)
            %������ "�������" ����� ����           
            semilogx((BUFTAU(g)+BUFTAU(g+1))/2,(BUFX(g)+BUFX(g+1))/2, ...
                'Marker','s','MarkerSize', 11,'MarkerEdgeColor','c', ...
                                              'MarkerFaceColor','c');
            hold on;
            kolk=kolk+1;%����������� ������� �����
        end
        if (imag(BUFLAM1(g))>0 && imag(BUFLAM1(g+1))<0) || ...
                (imag(BUFLAM2(g))>0 && imag(BUFLAM2(g+1))<0)
            %���������� � ������� ����������� ���������� ���������-����� (�������� � + �� -)
            %������ "�������" ����� ����            
            semilogx((BUFTAU(g)+BUFTAU(g+1))/2,(BUFX(g)+BUFX(g+1))/2, ...
                'Marker','s','MarkerSize', 11,'MarkerEdgeColor','c', ...
                                              'MarkerFaceColor','c');
            hold on;  
            kolk=kolk+1;%����������� ������� �����
%--------------------------------------------------------------------------            
        end                
    end
%������� grid on ������� ������������ ����� �� ������� ���    
grid on;
%��������� ����� ��� �������� ��� ������ title
title(sprintf('������ ����������� x%i (tau) ��� �=%i:',strx,at),'fontsize',15);
%�������� ����� �� ��� x � tau
xlabel('tau','fontsize',15); ylabel(sprintf('x%i ',strx),'fontsize',15,'rotation',0);
%������� ����������� ����� ������ �����    
fprintf('\n\t��� a=%i',at);
fprintf('\n\t����������� ������������ ���������� ����� ��� x%i(tau):%i',strx,stac);
fprintf('\n\t����������� ������������ ������������ ����� ��� x%i(tau):%i',strx,nstac);
fprintf('\n\t����������� ����� ������������ ���������� ��� x%i(tau):%i',strx,kolv);
fprintf('\n\t����������� ����� ���������� ���������-����� ��� x%i(tau):%i',strx,kolk);
fprintf('\n\t����������� ������ ����� ��� x%i(tau):%i',strx,ost);
fprintf('\n');
fig=fig+1;%����������� ������� fig
end