function [ecg1,i,energia] = criaecg( quant,ecg,pos,laco ) 
    ecg1=zeros(1,quant*360);
    for i=pos:(pos-1)+quant*360
        if laco ~= 1
            ecg1(i-(laco-1)*360)=ecg(i);
        end
        if laco == 1
            ecg1(i)=ecg(i);
        end
    end
    Ts=1/360;
    L=360*quant;
    To=L*Ts;
    t=(pos-1)*Ts:Ts:Ts*(laco*(L-1)+(laco-1));
    ecgfreq=fft(ecg1);
    k=0:L-1;
    o=k/To;
    energia=4.8*1.4*10^-3;
end
