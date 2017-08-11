clear all; close all; clc;
%addpath(genpath('../../'))
zoomquality = 'high'; %low or high
type = 'rectanglewedge'; %rectangle or rectanglewedge
gam_eV = 0.5e-3; 
%% LOAD
loadname = ['dipExciteOneWay_gam' strrep(num2str(gam_eV*1e3),'.','p') 'meV_B8T_' type '_' zoomquality 'reszoom.mat']; 
load(['../../output/dipole/' loadname])
%% PLOT SETTINGS
[cols,nams] = flatcolors;
coloutline = cols{8}*.6 + cols{4}*.4;
lwoutline = .5;

%% SORT OUTLINE
remcnct = cnct; pp = 1; newcnct{pp} = remcnct(1,:); remcnct(1,:) = []; 
while numel(remcnct)>1
	[row,col] = find(remcnct == newcnct{pp}(end));
    if ~isempty(row)
        newcnct{pp}(end+1) = remcnct(row,mod(col,2)+1); 
        remcnct(row,:) = []; 
    else
        pp = pp + 1; 
        newcnct{pp} = remcnct(1,:); 
        remcnct(1,:) = []; 
    end
end
%% PLOT OUTLINE

%Plot outline efficiently
for pp = 1:numel(newcnct)
    plot(nodes(newcnct{pp},1)*L*1e6, nodes(newcnct{pp},2)*L*1e6,'-b'); hold on 
    drawnow; pause(1)
end

