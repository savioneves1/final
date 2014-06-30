clear all
close all
clc

n=3;
M=textread('correntesphoto.txt','%f');
M=reshape(M,n,[])';
M(:,2)=-M(:,2);
figure(1);
tempo=0:0.00659631:length(M(:,2))*0.00659631;
tempo(length(tempo))=[];
plot(tempo,M(:,2),'b-o')
hold on;
M=textread('correntesdctphoto.txt','%f');
M=reshape(M,n,[])';
M(:,2)=-M(:,2);
plot(tempo,M(:,2),'r-x');
xlabel('Segundo (s)','FontSize',14.5);
ylabel('Corrente (A)','FontSize',14.5);
h_legend=legend('Normal','DCT');
set(h_legend,'FontSize',14.5);
grid on;