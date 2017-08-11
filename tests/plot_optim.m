clear all; close all; clc; 
[cols,nams] = flatcolors();
pars = {[1.001,inf,0,1],...
        [25,inf,0,1],...
        [8.7400,0.3646,-0.8589,2.2944],...
        [38.8435,0.5164,-4.8324,3.2006],...
        [53.7888,0.3964,-3.0917,3.6358]};
N = [150,350,190,250,250]; 
for pp = 1:numel(pars)
    [mesh{pp},nodes{pp}] = geomPeanut(pars{pp}(1),pars{pp}(2),pars{pp}(3),pars{pp}(4),N(pp),[],[],1);
    oscstrength_norm{pp} = optim_oscstrength(mesh{pp},1,[],0);
end

%% 
for pp = 1:numel(pars)
subaxis(numel(pars),2,1,pp)
    trimesh(mesh{pp}.t,mesh{pp}.p(:,1),mesh{pp}.p(:,2),'Color',cols{pp},'LineWidth',.2); %Plot with an off-black color [from flatcolors()]
    %Clean up
    axis equal off
    xlim(minmax(mesh{pp}.p(:,1))+[-1,1]*max(abs(mesh{pp}.p(:)))/15)
    ylim(minmax(mesh{pp}.p(:,2))+[-1,1]*max(abs(mesh{pp}.p(:)))/15)
    set(gca,'Fontsize',8,'LineWidth',.2); drawnow;
    title(['\alpha_n/||\Omega|| = ' num2str(max(oscstrength_norm{pp})*100,'%.1f') '%'])
    
subaxis(numel(pars),2,2,pp)
    fill(nodes{pp}(:,1),nodes{pp}(:,2),cols{pp})
    %Clean up
    axis equal off
    xlim(minmax(mesh{pp}.p(:,1))+[-1,1]*max(abs(mesh{pp}.p(:)))/15)
    ylim(minmax(mesh{pp}.p(:,2))+[-1,1]*max(abs(mesh{pp}.p(:)))/15)
    set(gca,'Fontsize',8,'LineWidth',.2); drawnow;
    title(['\alpha_n||\Omega|| = ' num2str(max(oscstrength_norm{pp})*100,'%.1f') '%, \zeta_n/||\Omega||^{1/2} = '])
end
