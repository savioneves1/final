%Camada PHY utiliza 6 bytes
%Camada MAC utiliza 25+txfifo bytes; txfifo pode possuir até 128 bytes; 

%System Specifications

%Radio
%Trasmit Mode
%Output power(dBm)   Current consumption(mA)
% 0                  17,4
%-1                  16,5
%-3                  15,2
%-5                  13,9
%-7                  12,5
%-10                 11,2
%-15                 9,9
%-25                 8,5
%Receive Mode
%Current Consumption(mA)=18,8
%Idle mode(microA)=426

function [bits, energ_trans ] = txpacote( ondauni )
    import java.lang.Integer java.lang.Float;
    
    txfifo=(128-6-25)*8;
    bits=0;
    for i=1:length(ondauni)
        bits=bits+length(Integer.toBinaryString(Float.floatToIntBits(ondauni(i))));
    end
    transmit=ceil(bits/txfifo);
    tempo_transmit=((transmit*((6+25)*8)+bits)/250000);
    energ_trans=tempo_transmit*17.4*(10^-3)*1.8 + (1-tempo_transmit)*1.8*426*(10^-6);
end