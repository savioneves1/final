clear all;
close all;
clc;

import java.lang.Integer java.lang.Float;

quant=1;
load 103_ch1_000to300.mat;
tamanho=size(ecg);
pos=1;
gerou_matuni=0;
gerou_matnouni=0;
energia_totuni=0;
energia_totnonuni=0;
energia_totnocr=0;

W=50;
q=0.5;
prob1=0.3;
prob2=0.4;
oldcenter=0;

energia_tot_transuni=0;
energia_tot_procuni=0;

energia_tot_transnocr=0;

energia_tot_transnonuni=0;
energia_tot_procnonuni=0;

energia_tot_transdct=0;
energia_tot_procdct=0;

energia_tot_sens=0;

bits_tot_nocr=0;
bits_tot_uni=0;
bits_tot_nonuni=0;
bits_tot_dct=0;

bits_nocr=0;
bits_uni=0;
bits_nonuni=0;
bits_dct=0;

coef_dct=zeros(1,360);

for i=1:round(tamanho(1)/360)
    i;
    %Sensing
    [ecg1,pos,energia_sens]=criaecg(quant,ecg,pos,i);
    pos=pos+1;
    
    %Processamento Uniforme
    if gerou_matuni==0
        [matbin]=gera_matbin(0.8);
        gerou_matuni=1;
    end
    if gerou_matuni~=0
        [ondauni,energia_procuni]=comp_sens_uni(ecg1,matbin);
        ondauni=ondauni';
    end
    

    %Processamento Não Uniforme
    if gerou_matnouni==0
        [matbinnon]=gera_matbinnon(0.8,W,q,prob1,prob2,ecg1);
        gerou_matnouni=1;
    end
    if gerou_matnouni~=0
        [maximo,minimo,center]=ROI(ecg1,W);
        if i==1
            oldcenter=center;
        end
        [ondanonuni,energia_procnonuni]=comp_sens_nonuni(ecg1,matbinnon,center,i,oldcenter);
        ondanonuni=ondanonuni';
        oldcenter=center;
    end
    
    
    
    %DCT
    cdct=dct(ecg1);
    energia_procdct=(0.05e-3*2.7^2)+(2.7*(1.196e-3*exp(2.7/(21.26*0.2)))*(0.08e6/191.42e6));
    if i==1
        coef_dct=cdct;
    else
        coef_dct=[coef_dct;cdct];
    end
     
    %Transmissão
    [bits_nocr,energia_transnocr]=txpacote(ecg1);
    [bits_uni,energia_transuni]=txpacote(ondauni);
    [bits_nonuni,energia_transnonuni]=txpacote(ondanonuni);
    
    %Energia total
    energia_totuni = energia_totuni + energia_transuni + energia_procuni + energia_sens;
    energia_totnocr = energia_totnocr + energia_sens + energia_transnocr;
    energia_totnonuni = energia_totnonuni + energia_transnonuni + energia_procnonuni + energia_sens;

    %bits
    bits_tot_nocr=bits_tot_nocr+bits_nocr;
    bits_tot_uni=bits_tot_uni+bits_uni;
    bits_tot_nonuni=bits_tot_nonuni+bits_nonuni;
    
    %energia compressão uniforme
    energia_tot_transuni=energia_tot_transuni+energia_transuni;
    energia_tot_procuni=energia_tot_procuni+energia_procuni;
    
    %energia sem compressão
    energia_tot_transnocr=energia_tot_transnocr+energia_transnocr;

    %energia compressão não uniforme
    energia_tot_transnonuni=energia_tot_transnonuni+energia_transnonuni;
    energia_tot_procnonuni=energia_tot_procnonuni+energia_procnonuni;

    %energia compressão dct
    energia_tot_procdct=energia_tot_procdct+energia_procdct;
    
    %energia sensing
    energia_tot_sens=energia_tot_sens+energia_sens;
end

dim=size(coef_dct);
varia=zeros(1,dim(2));
for i=1:dim(2)
    varia(1,i)=var(coef_dct(:,i));
end
[varia,index]=sort(varia,'descend'); 
energia_var=((0.67e-3*2.7^2)+2.7*(1.196e-3*exp(2.7/(21.26*0.2)))*(0.97e6/191.42e6))*360;
ecgdct=zeros(dim(1),length(varia)/2);
for i=1:length(varia)/2
    ecgdct(:,i)=coef_dct(:,index(i));
end
dim=size(ecgdct);
for i=1:dim(1)
    [bits_dct,energia_transdct]=txpacote(ecgdct(i,:));
    bits_tot_dct=bits_tot_dct+bits_dct;
    energia_tot_transdct=energia_tot_transdct+energia_transdct;
end
energia_totdct= energia_tot_transdct + energia_tot_procdct + energia_tot_sens+energia_var;

figure(1);
metodo=[1 2 3 4];
energia1=[energia_totnocr 0 0  0];
energia2=[0 energia_totuni 0  0];
energia3=[0 0 energia_totnonuni 0];
energia6=[0 0 0 energia_totdct];
bar(metodo,energia1,0.5,'b');
hold on;
bar(metodo,energia2,0.5,'r');
bar(metodo,energia3,0.5,'g');
bar(metodo,energia6,0.5,'c');
hold off;
h_legend=legend('Sem compressão','Compressed Sensing Uniform','Compressed Sensing Non-Uniform','DCT');
set(h_legend,'FontSize',14.5);
xlabel('Método','FontSize',14.5);
ylabel('Consumo de Energia(J)','FontSize',14.5);
grid on;

figure(2);
pie([energia_tot_transnocr,energia_tot_sens]);
h_legend4=legend('Transmissão','Aquisição');
set(h_legend4,'FontSize',14.5);

figure(3);
pie([energia_tot_transuni,energia_tot_procuni,energia_tot_sens]);
h_legend5=legend('Transmissão','Processamento','Aquisição');
set(h_legend5,'FontSize',14.5);

figure(4);
pie([energia_tot_transnonuni,energia_tot_procnonuni,energia_tot_sens]);
h_legend6=legend('Transmissão','Processamento','Aquisição');
set(h_legend6,'FontSize',14.5);

figure(5);
pie([energia_tot_transdct,energia_tot_procdct,energia_var,energia_tot_sens]);
h_legend7=legend('Transmissão','Processamento-DCT','Processamento-Variância','Aquisição');
set(h_legend7,'FontSize',14.5);

