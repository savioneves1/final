clear all;
close all;
clc;

n=12;%Número de colunas no arquivo de temperatura
COL=4;%Número de colunas a serem lidas no arquivo de temperatura
n2=14;%Número de colunas no arquivo de coeficientes
COL2=6;%Número de colunas no arquivo de coeficientes

M=textread('dadosdctphoto_photo.txt','%2c');
M=hex2dec(char(M));
M=reshape(M,n,[])';
dim=size(M);

M2=textread('dadosdctphoto_coef.txt','%2c');
M2=hex2dec(char(M2));
M2=reshape(M2,n2,[])';
dim2=size(M2);

elem=zeros(dim(1),COL);
for i=1:COL
    elem(:,i)=M(:,8+i);
end
elem2=zeros(dim2(1),COL2);
for i=1:COL2
    elem2(:,i)=M2(:,8+i);
end
for i=1:2:COL
    elem(:,i)=bitshift(elem(:,i),8);
end
for i=1:2:COL2
    elem2(:,i)=bitshift(elem2(:,i),8);
end
meio=zeros(dim(1),COL/2);
for i=1:COL/2
    for j=1:dim(1)
        meio(j,i)=elem(j,2*i)+elem(j,2*i-1);
    end
end
meio2=zeros(dim2(1),COL2/2);
for i=1:COL2/2
    for j=1:dim2(1)
        meio2(j,i)=elem2(j,2*i)+elem2(j,2*i-1);
    end
end

fim=meio;
fim(:,1)=typecast(uint16(fim(:,1)),'int16');
dim=size(fim);
grupo_maior=fim(dim(1),2);
bits=zeros(grupo_maior,1);
bits(:,1)=14*8;
grupo_atual=1;
j=1;
for i=1:grupo_maior
    count = 0;
    while ((j<=dim(1))&&(fim(j,2) == grupo_atual))
        count=count+1;
        j=j+1;
    end
    bits(grupo_atual,1)=bits(grupo_atual,1)+count*12*8;
    grupo_atual=grupo_atual+1;
end

bits_normal=10*8*32;
metodo=[1 2];
bits_normal=[bits_normal 0];
media_fdct=[0 mean(bits)];
media=[bits_normal(1) mean(bits)];
err=[0 (std(bits)/sqrt(10))*2.262];
figure(1);
bar(metodo,bits_normal,0.5,'b');
hold on;
bar(metodo,media_fdct,0.5,'r');
errorbar(metodo,media,err,'kx');
hold off;
xlabel('Método','FontSize',11.5);
ylabel('Bits','FontSize',11.5);
h_legend=legend('Normal','DCT');
set(h_legend,'FontSize',11.5);
bar_bits_photo = applyhatch_pluscolor(gcf,'|-',0,0);
imwrite(bar_bits_photo,'bar_bits_photo.png','png');