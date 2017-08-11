omegac = e*model.B*vf^2/(model.ef_eV*ev2jo); %Classical cyclotron: graphene
sigmac = conducIntra(omegac*hbar_eV,model.ef_eV); %Drude conductivity evaluated at cyclotron freq.
zetac = 2i*eps0*omegac*model.L/sigmac;


kk=1;
k = kvec(kk,:);
[~,eW] = evalc('calcLatticeCoulombOptimized(blochmesh,k)'); %Run without textual output
[~,eDRi,eDL] = evalc('calcLatticeDifferentialIsotropic(blochmesh,k)'); %Run without textual output
[~,eDRa] = evalc('calcLatticeDifferentialAnisotropic(blochmesh,k)'); %Run without textual output

%%



set_figsize([],50,25);
boti = find(eneeig_eV{kk}/eneB_eV > 0 & eneeig_eV{kk}/eneB_eV < 1); topi = find(eneeig_eV{kk}/eneB_eV > 1); 
eelist = [boti(end-20:end); topi(1:(36-numel(boti)))];
for eee=1:36
    ee = eelist(eee);
    
    subaxis(4,9,eee,'M',0.01,'MT',.05,'SH',.01,'SV',.05)

    lambda = eneeig_eV{kk}/eneB_eV; 
    A = (1i/(2*pi*zetac))*(eDL\(eDRa*eW)) - lambda(ee)*((1/(2*pi*zetac))*(eDL\(eDRi*eW)) + eye(size(eW))) + lambda(ee)^3*eye(size(eW));
    rho = eW*null(A);
    trisurf(blochmesh.remesh.t,blochmesh.remesh.p(:,1),blochmesh.remesh.p(:,2),real(blochmesh.remesh.values(rho)),'EdgeColor','none');
    hold on
    plot([-.75,.25,.75,-.25,-.75],[-1,-1,1,1,-1]*R(2,2)/2,'--','color',cols{3}); 
    hold off
    
    cvals=sort(abs(real(rho(:))),'descend');
    caxis([-1,1]*cvals(20))
    %shading interp;
    colormap('bluewhitered_mod'); 
    freezeColors;
    view(2);
    axis equal off;
    xlim([-1,1]*1.35); ylim([-1,1]*1.2)
    title(['\omega/\omega_c = ' num2str(real(eneeig_eV{kk}(ee)/eneB_eV),2)]);
    drawnow;% pause(.25)
end