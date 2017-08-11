clc; clear all; close all; 

addpath(genpath('..'))
%Mesh settings

Ne = 60;
Ncirc0 = 60;

latticetype = 'square';
meshfun = 'const';
adivdloop = 1.5;[6,5,4:-.5:1.5];
for adivd = adivdloop
    %Inclusion boundary
    Ncirc = round(Ncirc0*2/adivd);
    theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1);
    inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';
    
    %Meshing
    blochmesh = geomPeriodicInclusion(latticetype,inclusion,Ne,meshfun);
end


