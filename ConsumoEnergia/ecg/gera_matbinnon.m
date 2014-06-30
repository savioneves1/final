function [ matbin ] = gera_matbinnon( CR,W,q,prob1,prob2,onda )
    
    matbin1=binornd(1,prob1,q*360*(1-CR),W+1);
    matbin2=binornd(1,prob2,360*(1-CR)-(q*360*(1-CR)),360);
    
    [maximo,minimo,center]=ROI(onda,W);
    MATRIZ=zeros(q*360*(1-CR),360);
    for i=minimo:maximo
        for j=1:q*360*(1-CR)
            MATRIZ(j,i)=matbin1(j,i-(minimo-1));
        end
    end
    matbin=[MATRIZ ; matbin2];
end