clc,clear,close all
set(0,'DefaultFigureWindowStyle','docked')
ps=readmatrix("Treloar_pureshear.csv");

kgfcm2ToNmm2 = 9.80665/(10)^2;
x = ps(:,1);
y = ps(:,2)*kgfcm2ToNmm2;
plot(x,y,'ko-')
grid on