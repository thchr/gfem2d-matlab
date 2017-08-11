clear all; close all; clc;
set_figsize(1,35,30)
B = [0,4,8];
for kk = [3,2,1,4,5,6]
    
    load(['../output/ribbons/scissored_july/scmag_triangular_numcell20_cut40_Ns70_Ncirc36_savephi0_' num2str(kk) '.mat'])
    cols = flatcolors; models(1).B = 0;
    
 %   if kk == 3;
        for mm = 2:3
            etemp = real(eneeig_eV{mm}(:,1));
            kill{mm} = etemp(etemp<=omegac_eV(B(mm)) & etemp > 1e-10);
            nzeros(mm) = nnz(etemp < 1e-10);
            eneeig_eV{mm}(1:nzeros(mm),:) = NaN;
        end
  %  end
    
    tol = .5e-4;
    for mm = 2:3
        for kkk = 1:size(kvec,1)
            for ee = nzeros(mm):nzeros(mm)+round(numel(kill{mm})*1.2)
                for eek = 1:numel(kill{mm})
                    if real(eneeig_eV{mm}(ee,kkk)) >= kill{mm}(eek)-tol && real(eneeig_eV{mm}(ee,kkk)) <= kill{mm}(eek)+tol
                        eneeig_eV{mm}(ee,kkk) = NaN;
                    end
                end
            end
        end
    end
    
    for mm = 1:3
        subaxis(1,3,mm)
        plot(kvec(:,1),real(eneeig_eV{mm}),'.k'); hold on
        plot([0,pi],[1,1]*omegac_eV(models(mm).B),'-','color',cols{8},'LineWidth',2); hold off
        view(2); ylim([0,.12]); xlim([0,pi]); set(gca,'Layer','Top'); box on
        drawnow
    end
end

