function M=Main
%Очищаем Окно Команд и область переменных
clc
clear all;
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

x1j=x1;%Присваиваем значение x1 в буферную переменную
x2j=x2;%Присваиваем значение x2 в буферную переменную
%Откроем файл OUTValues.dat для дальнейшей записи таблицы данных в него
outfile=fopen('D:\Мат.модели\Курсач\OUTValues.dat','w');
%--------------------------------------------------------------------------
%Блок аналитических преобразований
fprintf('Умножим (1) уравнение на p4:');
temp1=dx1dt*p4
fprintf('Вычтем (1) уравнение из (2):');
temp2=dx2dt-temp1
fprintf('Выразим x1:');
x1=solve(temp2,x1)
fprintf('Подставим x1 в (1):');
temp3=dx1dt;
temp3=subs(temp3,x1)
fprintf('Решим квадратное уравнение по переменной tau:');
taut=solve(temp3,tau)
fprintf('Строим матрицу Якоби:');
maj=jacobian([dx1dt,dx2dt],[x1j,x2j])
%--------------------------------------------------------------------------
%Основная часть (логика,цикл, вывод результатов и графиков)
%Создаем буферный массив 1*2 для записи корней уравнения (tau) в него
kor=[0 0];
%Создаем буферный массив 2*2 для записи матрицы Якоби в него
majt=[0 0;0 0];
%Создаем пустые массивы для вывода всех данных в виде таблицы
KP=[];A=[];X2=[];TAU=[];K=[];X1=[];LAM1=[];LAM2=[];
%Создаем пустые буферные массивы для вывода графиков по точкам
BUFA=[];BUFTAU=[];BUFX1=[];BUFX2=[];BUFLAM1=[];BUFLAM2=[];
GRAFA=[];GRAFTAU=[];GRAFX1=[];GRAFX2=[];GRAFLAM1=[];GRAFLAM2=[];
%Создаем пустые буферные массивы для сортировки точек графиков
BUFTAU1=[];BUFTAU2=[];BUFX11=[];BUFX12=[];BUFX21=[];BUFX22=[];
TEMPLAM11=[];TEMPLAM12=[];TEMPLAM21=[];TEMPLAM22=[];BUFA1=[];BUFA2=[];
iter=1;%Счетчик общего колличества точек
fig=1;%Переменная, отвечающая за номер фигуры выведенного графика
at=1;%Эквивалент переменной а, меняющий значения в цикле а=1 2 4
%Начало цикла, пока а не станет равной 4
while at<=4
    for x2t=-5.0:0.1:10.0
        try
        kor=double(subs(taut,[x2,a],[x2t,at]));
        catch
        %Блок try-catch необходим для фильтрации ошибок и предупреждений
        %вида 'Деление на ноль' при решении квадратного уравнения по tau
        end
        %Далее проверяем является ли корень tau комплексным числом
        if imag(kor)~=0
            kor=abs(kor);
        end
        %Пока выполняется цикл будем выводить результаты в Окно Команд 
        fprintf('a=%d x2=%4.1f\n',at,x2t);
        %Так как у нас мб 2 значения tau=>цикл по kor(k) (или tau)
        for k=1:length(kor) %Длина kor принимает значения от 1 до 2
            %Подставляем значения a,x2,tau в буферную переменную x1t
            x1t=double(subs(x1,[a,x2,tau],[at,x2t,kor(k)]));
            %Подставляем значения a,x1,x2,tau в буферную матрицу Якоби majt
            majt=double(subs(maj,[a,x1j,x2j,tau],[at,x1t,x2t,kor(k)]));
            %Вычисляем сумму элементов главной диагонали матрицы majt
            sd=majt(1,1)+majt(2,2);
            %Находим вектор собственных значений матрицы Якоби majt
            ls=double(eig(majt));
            %Для вывода данных в общую таблицу постоянно записываем все 
            %значения в массивы 
            A(end+1)=at;
            X2(end+1)=x2t;
            KP(end+1)=iter;
            %Увеличиваем счетчик общего колличества точек
            iter=iter+1;
            %Вводим условие из задания
            if (kor(k)>0 && x1t>0 && x1t<1)
            %Пока выполняется цикл будем выводить результаты в Окно Команд 
                 fprintf('\ttau(%d)= %9.5f | x1=%9.5f | L1=%s | L2=%s \n', ...
                 k,kor(k),x1t,num2str(ls(1)),num2str(ls(2)));
            %Далее логически разделяем поэлементно значения в массивах
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
            %Если же корни tau отсутствуют, присваиваем пустое значение
            else
                 TAU(end+1)=NaN;
                 K(end+1)=NaN;
                 X1(end+1)=NaN;
                 LAM1(end+1)=NaN;
                 LAM2(end+1)=NaN;                 
                 fprintf('\t%d-й корень tau отсутствует\n',k);
            continue
            %Если k-ый корень отсутствует, понижаем счетчик кол-ва точек
                 iter=iter-1;
            end
        end
        fprintf('\n');
    end 
    %Дабы избежать "зиг-загообразности" точек на графиках:
    %Переворачиваем ранее записанный массив вторых корней уравнения
    %и сопутствующие ему буферные массивы значений
    BUFA2=fliplr(BUFA2);
    BUFTAU2=fliplr(BUFTAU2);
    BUFX12=fliplr(BUFX12);
    BUFX22=fliplr(BUFX22);
    TEMPLAM21=fliplr(TEMPLAM21);
    TEMPLAM22=fliplr(TEMPLAM22);
    %И скрепляем с массивами первых корней уравнения
    BUFA=[BUFA1,BUFA2];
    BUFTAU=[BUFTAU1,BUFTAU2];
    BUFX1=[BUFX11,BUFX12];
    BUFX2=[BUFX21,BUFX22];
    BUFLAM1=[TEMPLAM11,TEMPLAM21];
    BUFLAM2=[TEMPLAM12,TEMPLAM22];
    %Передаем значения другим массивам для вывода точек всех графиков
    GRAFA=[GRAFA,BUFA];
    GRAFTAU=[GRAFTAU,BUFTAU];
    GRAFX1=[GRAFX1,BUFX1];
    GRAFX2=[GRAFX2,BUFX2];
    GRAFLAM1=[GRAFLAM1,BUFLAM1];
    GRAFLAM2=[GRAFLAM2,BUFLAM2];
    %Используем собственные функции построения графиков для x1 и x2
    fig=mylogplot(at,fig,BUFTAU,BUFLAM1,BUFLAM2,BUFX1);
    fig=mylogplot(at,fig,BUFTAU,BUFLAM1,BUFLAM2,BUFX2);
    %Обнуляем массивы, чтобы точки для разных значений а
    %не лежали в одной плоскости
    BUFA=[];BUFTAU=[];BUFX1=[];BUFX2=[];BUFLAM1=[];BUFLAM2=[]; 
    BUFTAU1=[];BUFTAU2=[];BUFX11=[];BUFX12=[];BUFX21=[];BUFX22=[];
    TEMPLAM11=[];TEMPLAM12=[];TEMPLAM21=[];TEMPLAM22=[];BUFA1=[];BUFA2=[];
    %Умножаем фактическое а на 2
     at=at*2;
end
%Инициализируем таблицу значений всех точек и выводим в Окно Команд 
T=table;
T.N=KP';T.a=A';T.x2=X2';T.tau=TAU';
T.k=K';T.x1=X1';T.LAM1=LAM1';T.LAM2=LAM2'
%Инициализируем таблицу значений точек по которым строим графики
G=table;
G.N=(1:length(GRAFTAU))';G.a=GRAFA';G.x2=GRAFX2';G.tau=GRAFTAU';
G.x1=GRAFX1';G.LAM1=(double(GRAFLAM1))';G.LAM2=(single(GRAFLAM2))'
%Строим в отдельном окне легенду точек, общую для всех графиков: 
figure(fig);
imshow('D:\Мат.модели\Курсач\Legend.png');
%Записываем таблицу в отдельный файл OUTValues.dat 
writetable(T,'D:\Мат.модели\Курсач\OUTValues.dat');
%Закрываем файл
fclose(outfile);
%Блок проверки
Test(GRAFA,GRAFTAU,GRAFX1,GRAFX2);
end