function [fig]=mylogplot(at,fig,BUFTAU,BUFLAM1,BUFLAM2,BUFX)
kolv=0;%Счетчик точек вещественной бифуркации
kolk=0;%Счетчик точек бифуркации Андронова-Хопфа
stac=0;%Счетчик стационарных устойчивых точек
nstac=0;%Счетчик стационарных неустойчивых точек
ost=0;%Счетчик точек, не входящих ни в одно условие
%В зависимости от входного массива задаем строковой переменной strx
%значение 1 или 2
if (rem(fig,2)==1)
    strx=1;
elseif (rem(fig,2)==0)
    strx=2;
end
    %Строим график в новом окне с индексом fig
    figure(fig);
    %Цикл с условием устранения выхода за пределы массива length(BUFTAU)-1
    for g=1:length(BUFTAU)-1
    %Вводим условия для 4 типов точек и обозначаем их на графиках
    %в логарифмическом масштабе по оси tau
        if (real(BUFLAM1(g))<0 && real(BUFLAM2(g))<0)
            %Стационарные устойчивые точки
            %Маркер "звездочка" зеленого цвета
            semilogx(BUFTAU(g),BUFX(g),'Marker','p','MarkerEdgeColor','g', ...
                                                    'MarkerFaceColor','g');
            %Команда hold on включает режим сохранения текущего графика 
            hold on;
            stac=stac+1;%Увеличиваем счетчик точек
%--------------------------------------------------------------------------            
        elseif ((real(BUFLAM1(g))>0 && real(BUFLAM2(g))>0) || ...
                (real(BUFLAM1(g))<0 && real(BUFLAM2(g))>0) || ...
                (real(BUFLAM1(g))>0 && real(BUFLAM2(g))<0))
            %Стационарные неустойчивые точки
            %Маркер "o" красного цвета            
            semilogx(BUFTAU(g),BUFX(g),'Marker','o','MarkerEdgeColor','r', ...
                                                    'MarkerFaceColor','r');
            hold on;    
            nstac=nstac+1;%Увеличиваем счетчик точек
%--------------------------------------------------------------------------            
        else
           %Другие точки
           %Маркер "о" черного цвета
           semilogx((BUFTAU(g))',(BUFX(g))','ko');
           hold on;
           ost=ost+1;%Увеличиваем счетчик точек
        end
%----------------------------------------------------------------------------------------------
		if (real(BUFLAM1(g))<0 && real(BUFLAM1(g+1))>0) || ...
                (real(BUFLAM2(g))<0 && real(BUFLAM2(g+1))>0) 
            %Промежуток с точками вещественной бифуркации (меняется с - на +)
            %Маркер "звезда" синего цвета
            semilogx((BUFTAU(g)+BUFTAU(g+1))/2,(BUFX(g)+BUFX(g+1))/2, ...
                'Marker','h','MarkerSize',11,'MarkerEdgeColor','b', ...
                                             'MarkerFaceColor','b');
            hold on;
            kolv=kolv+1;%Увеличиваем счетчик точек
        end
        if (real(BUFLAM1(g))>0 && real(BUFLAM1(g+1))<0) || ...
                (real(BUFLAM2(g))>0 && real(BUFLAM2(g+1))<0)
            %Промежуток с точками вещественной бифуркации (меняется с + на -)
            %Маркер "звезда" синего цвета            
            semilogx((BUFTAU(g)+BUFTAU(g+1))/2,(BUFX(g)+BUFX(g+1))/2, ...
                'Marker','h','MarkerSize',11,'MarkerEdgeColor','b', ...
                                             'MarkerFaceColor','b');
            hold on;    
            kolv=kolv+1;%Увеличиваем счетчик точек
        end
%--------------------------------------------------------------------------        
		if (imag(BUFLAM1(g))<0 && imag(BUFLAM1(g+1))>0) || ...
                (imag(BUFLAM2(g))<0 && imag(BUFLAM2(g+1))>0) 
            %Промежуток с точками комплексной бифуркации Андронова-Хопфа (меняется с - на +)
            %Маркер "квадрат" цвета циан           
            semilogx((BUFTAU(g)+BUFTAU(g+1))/2,(BUFX(g)+BUFX(g+1))/2, ...
                'Marker','s','MarkerSize', 11,'MarkerEdgeColor','c', ...
                                              'MarkerFaceColor','c');
            hold on;
            kolk=kolk+1;%Увеличиваем счетчик точек
        end
        if (imag(BUFLAM1(g))>0 && imag(BUFLAM1(g+1))<0) || ...
                (imag(BUFLAM2(g))>0 && imag(BUFLAM2(g+1))<0)
            %Промежуток с точками комплексной бифуркации Андронова-Хопфа (меняется с + на -)
            %Маркер "квадрат" цвета циан            
            semilogx((BUFTAU(g)+BUFTAU(g+1))/2,(BUFX(g)+BUFX(g+1))/2, ...
                'Marker','s','MarkerSize', 11,'MarkerEdgeColor','c', ...
                                              'MarkerFaceColor','c');
            hold on;  
            kolk=kolk+1;%Увеличиваем счетчик точек
%--------------------------------------------------------------------------            
        end                
    end
%Команда grid on наносит координатную сетку на текущие оси    
grid on;
%Размещаем текст над графиком при помощи title
title(sprintf('График зависимости x%i (tau) при а=%i:',strx,at),'fontsize',15);
%Помещаем текст на оси x и tau
xlabel('tau','fontsize',15); ylabel(sprintf('x%i ',strx),'fontsize',15,'rotation',0);
%Выводим колличество точек разных типов    
fprintf('\n\tПри a=%i',at);
fprintf('\n\tКолличество стационарных устойчивых точек для x%i(tau):%i',strx,stac);
fprintf('\n\tКолличество стационарных неустойчивых точек для x%i(tau):%i',strx,nstac);
fprintf('\n\tКолличество точек вещественной бифуркации для x%i(tau):%i',strx,kolv);
fprintf('\n\tКолличество точек бифуркации Андронова-Хопфа для x%i(tau):%i',strx,kolk);
fprintf('\n\tКолличество других точек для x%i(tau):%i',strx,ost);
fprintf('\n');
fig=fig+1;%Увеличиваем счетчик fig
end