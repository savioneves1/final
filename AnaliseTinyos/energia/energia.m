clear all
close all
clc

energ=zeros(1,10);
energ2=zeros(1,10);
tensao=2.64;
interval=0.00659631;
n=3;
M=textread('correntesphoto.txt','%f');
M=reshape(M,n,[])';
dim=size(M);
M(:,2)=-M(:,2);
M2=textread('correntesdctphoto.txt','%f');
M2=reshape(M2,n,[])';
dim=size(M2);
M2(:,2)=-M2(:,2);

tempo=0:interval:length(M(:,2))*interval;
tempo(length(tempo))=[];
index=1;
for i=1:10
    while (tempo(index)<=32*i)
        energ(i)=energ(i)+M(index,2)*tensao*interval;
        energ2(i)=energ2(i)+M2(index,2)*tensao*interval;
        if index==length(tempo)
            break;
        end
        if index<length(tempo)
            index=index+1;
        end
    end
end

metodo =[1 2];
media1=[mean(energ) 0];
media2=[0 mean(energ2)];
media=[mean(energ) mean(energ2)];
err=[(std(energ)/sqrt(10))*2.262 (std(energ2)/sqrt(10))*2.262];
figure(1);
bar(metodo, media1,0.5,'b');
hold on;
bar(metodo,media2,0.5,'r');
errorbar(metodo,media,err,'kx');
hold off;
xlabel('Método','FontSize',11.5);
ylabel('Energia por Ciclo de Transmissão(J)','FontSize',11.5);
h_legend=legend('Normal','DCT');
set(h_legend,'FontSize',11.5);
bar_ene_photo = applyhatch_pluscolor(gcf,'|-',0,0);
imwrite(bar_ene_photo,'bar_ene_photo.png','png');

energ=zeros(1,10);
energ2=zeros(1,10);
tensao=2.64;
interval=0.00659631;
n=3;
M=textread('correntestemp.txt','%f');
M=reshape(M,n,[])';
dim=size(M);
M(:,2)=-M(:,2);
M2=textread('correntesdcttemp.txt','%f');
M2=reshape(M2,n,[])';
dim=size(M2);
M2(:,2)=-M2(:,2);

tempo=0:interval:length(M(:,2))*interval;
tempo(length(tempo))=[];
index=1;
for i=1:10
    while (tempo(index)<=33*i)
        energ(i)=energ(i)+M(index,2)*tensao*interval;
        energ2(i)=energ2(i)+M2(index,2)*tensao*interval;
        if index==length(tempo)
            break;
        end
        if index<length(tempo)
            index=index+1;
        end
    end
end

metodo =[1 2];
media1=[mean(energ) 0];
media2=[0 mean(energ2)];
media=[mean(energ) mean(energ2)];
err=[(std(energ)/sqrt(10))*2.262 (std(energ2)/sqrt(10))*2.262];
figure(3);
bar(metodo, media1,0.5,'b');
hold on;
h_bar=bar(metodo,media2,0.5,'r');
errorbar(metodo,media,err,'kx');
hold off;
xlabel('Método','FontSize',11.5);
ylabel('Energia por Ciclo de Transmissão(J)','FontSize',11.5);
h_legend=legend('Normal','DCT');
set(h_legend,'FontSize',11.5);
bar_ene_temp = applyhatch_pluscolor(gcf,'|-',0,0);
imwrite(bar_ene_temp,'bar_ene_temp.png','png');