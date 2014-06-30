clear all
close all
clc

load('PHOTODCT.mat');
LUZdct=T;
load('LUZ.mat');
LUZnormal=luz';

dim1=size(LUZdct);
dim2=size(LUZnormal);
menor=min(dim1(2),dim2(2));

erro=zeros(1,menor);
for i=1:menor
    erro(i)=(LUZdct(i)-LUZnormal(i))^2;
end
erromed=mean(erro);

plot(1:menor,LUZdct(1:menor),'r-x');
hold on
plot(1:menor,LUZnormal(1:menor),'b-o');
hold off
xlabel('Segundo (s)','FontSize',14.5);
ylabel('Tensão (V)','FontSize',14.5);
h_legend=legend('DCT','Normal');
set(h_legend,'FontSize',14.5);
