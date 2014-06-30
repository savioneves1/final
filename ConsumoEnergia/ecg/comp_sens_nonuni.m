function [ ondacompr,energia_procnonuni] = comp_sens_nonuni( onda,matriz,center,laco,oldcenter )
    CR=0.8;
    q=0.5;
    tamanho=size(matriz);
    if laco == 1
        onda=onda';
        ondacompr=matriz*onda;
    end
    if laco ~= 1
        mat1=matriz(1:q*360*(1-CR),:);
        mat2=matriz(q*360*(1-CR)+1:tamanho(1),:);
        diferenca=center-oldcenter;
        mat3=circshift(mat1,[0,diferenca]);
        MATRIZ=[mat3;mat2];
        onda=onda';
        ondacompr=MATRIZ*onda;
    end
    energia_procnonuni=(0.67e-3*2.7^2)+(2.7*(1.196e-3*exp(2.7/(21.26*0.2)))*(0.97e6/191.42e6));
end