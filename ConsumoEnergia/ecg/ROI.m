function [ maximo,minimo,center ] = ROI( onda,W )
    [valor,center]=max(onda);
    maximo=center+(W/2);
    minimo=center-(W/2);
end

