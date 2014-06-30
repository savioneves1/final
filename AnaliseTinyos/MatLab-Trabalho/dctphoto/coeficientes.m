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

%Escrita no arquivo
fileID=fopen('dadosdctphoto_dct32.txt','w');
dim=size(fim);
grupo_atual=1;
coef_atual=meio2(1,1)+bitshift(meio2(1,2),16);
bit_atual=1;
grupo_maior=fim(dim(1),2);
fprintf(fileID,'%d\n',grupo_maior);

for i=1:dim(1)
    while((bitget(coef_atual,bit_atual)~=1)&&(bit_atual<=32))
        bit_atual=bit_atual +1;
        fprintf(fileID,'0 ');
    end
    
    if fim(i,2) == grupo_atual
        fprintf(fileID,'%d ',fim(i,1));
        bit_atual=bit_atual+1;
    end
    if fim(i,2) ~= grupo_atual
        grupo_atual=fim(i,2);
        coef_atual=meio2(grupo_atual,1)+bitshift(meio2(grupo_atual,2),16);
        bit_atual=2;
        fprintf(fileID,'\n%d ',fim(i,1));
    end
    
    if (i==dim(1)&&(bit_atual<=32))
        while(bit_atual<=32)
            bit_atual=bit_atual +1;
            fprintf(fileID,'0 ');
        end
    end
end
fclose(fileID);