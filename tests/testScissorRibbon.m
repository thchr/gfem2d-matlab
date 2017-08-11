clear all;close all; clc

addpath(genpath('..'))

latticetype = 'triangular'; %Lattice type (square is the alternative)

numcell = 21; %Number of repeated unit cells in ribbon (vertical dir.)

Ns = 53; %Number of points along each unit cell side

Ncirc = 36; %Number of points used in circular inclusion

adivd = 2; %Period divided by inclusion diameter


theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----
yincl = sort(inclusion(:,2)); ystep = yincl(yincl > 0.001)/sind(60); ystep = ystep(1:2:15);
cut = .5 - ystep(2); %Reduction of edge unit cells in percent along R2

mesh = geomRibbonInclusionScissor(latticetype, numcell,inclusion,cut,Ns,...
    @(x,y) 10*.015*250/Ns*ones(size(x)),1);
hold on
plot(inclusion(:,1)-(numcell-1)/2*cosd(60),inclusion(:,2)-(numcell-1)/2*sind(60),'-r.')
hold off
%xlim([-.6,.6]-(numcell-1)/2*cosd(60));ylim([-.6,.6]-(numcell-1)/2*sind(60));
%set(get(handle(gcf),'JavaFrame'),'Maximized',1);
