clear all
close all
clc

load('DCTTEMPERATURA.mat');
Tdct=T;
load('temperatura.mat');
Tnormal=T;

dim1=size(Tdct);
dim2=size(Tnormal);
menor=min(dim1(1),dim2(1));

erro=zeros(1,menor);
for i=1:menor
    erro(i)=(Tdct(i)-Tnormal(i))^2;
end
erromed=mean(erro);

plot(1:menor,Tdct(1:menor),'r-x');
hold on
plot(1:menor,Tnormal(1:menor),'b-o');
hold off

xlabel('Segundo (s)','FontSize',14.5);
ylabel('Temperatura (ºC)','FontSize',14.5);
h_legend=legend('DCT','Normal');
set(h_legend,'FontSize',14.5);