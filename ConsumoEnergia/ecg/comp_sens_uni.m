function [ ondacompr,energia ] = comp_sens_uni( onda,matbin )
    onda=onda';
    ondacompr=matbin*onda;
    energia=(0.67e-3*2.7^2)+(2.7*(1.196e-3*exp(2.7/(21.26*0.2)))*(0.97e6/191.42e6));
end