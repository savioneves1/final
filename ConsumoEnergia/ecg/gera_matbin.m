function [ matbin ] = gera_matbin( CR )
    prob=0.3;
    matbin=binornd(1,prob,360*(1-CR),360);
end