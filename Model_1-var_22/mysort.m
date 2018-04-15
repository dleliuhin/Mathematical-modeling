function [BUFTAU,BUFX1,BUFX2,BUFLAM1,BUFLAM2]=mysort(BUFTAU,BUFX1,BUFX2,BUFLAM1,BUFLAM2)
for m=1:(length(BUFTAU)-1)
    for q=1:(length(BUFTAU)-m)
        if (BUFTAU(q)>BUFTAU(q+1))
            buftau=BUFTAU(q);
            bufx1=BUFX1(q);
            bufx2=BUFX2(q);
            buflam1=BUFLAM1(q);
            buflam2=BUFLAM2(q);
            
            BUFTAU(q)=BUFTAU(q+1);
            BUFX1(q)=BUFX1(q+1);
            BUFX2(q)=BUFX2(q+1);
            BUFLAM1(q)=BUFLAM1(q+1);
            BUFLAM2(q)=BUFLAM2(q+1);
            
            BUFTAU(q+1)=buftau;
            BUFX1(q+1)=bufx1;
            BUFX2(q+1)=bufx2;
            BUFLAM1(q+1)=buflam1;
            BUFLAM2(q+1)=buflam2;
        end
    end
end

end