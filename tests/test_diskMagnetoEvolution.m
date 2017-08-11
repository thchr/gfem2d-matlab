clc; clear all; close all;
ConstantsUnits0; cols = flatcolors(); Lcols = cols([14,21,3]);

%% Mesh
mesh = geomEllipse([1,1],175,.045,1);
L = [100,200,400]*1e-9;

%% System matrices
[mesh.DRi,mesh.DL] = calcDifferential(mesh.p,mesh.t);
mesh.DRa = calcDifferential(mesh.p,mesh.t,[0,1;-1,0]);
mesh.V = calcCoulomb(mesh.p,mesh.t);

%% Graphene
ef = .2;
gam = 1e-3;
ene = linspace(0.001,.1,100);
B = linspace(0.01,7,10);
omegac = e*B*vf^2/(ef*ev2jo)*hbar_eV;

%% Model setup
models.type = 'magneto';
models.system = 'graphene';
models.ef_eV = ef;
Emax = 2;
for ll = 1:numel(L)
    models.L = L(ll);
    
    for bb=1:numel(B)
        models.B  = B(bb);
        models.gamc = gam/omegac(bb);
        models.keep_eV = [omegac(bb)+1e-4,ef];
        eigs{bb,ll} = calcEigenEnergiesAny(mesh,models);
        for ee = 1:Emax
            sorteigs(ee,bb,ll) = eigs{bb,ll}(ee);
        end
    end
    
    %% Boring setup (no magnetic field)
    eigs0{ll} = calcEigenEnergiesAny(mesh,struct('type','isotropic','approx','intra','ef_eV',ef,'L',L(ll),'gam_eV',gam,'keep_eV',[0,ef]));
    eigs0_noloss{ll} = calcEigenEnergiesAny(mesh,struct('type','isotropic','approx','intra','ef_eV',ef,'L',L(ll),'keep_eV',[0,ef]));
    eigs0{ll} = eigs0{ll}(1:Emax); eigs0_noloss{ll} = eigs0_noloss{ll}(1:Emax);
for ll=1:numel(L)
    %% Guess at perturbation
    prefactor = .475;
    delta{ll} = bsxfun(@times,[-1;1].*(1-1i*gam./eigs0_noloss{ll}/2),omegac*prefactor);
    delta_full{ll} = bsxfun(@times,[-1;1].*eigs0{ll}./(2*eigs0{ll}+1i*gam)*2,omegac*prefactor);
end

%% Plot

ylimsRe = [0,.1];
xlims = minmax(B)+[-1,1]*.25; 
markerstyle = {'o','MarkerFaceColor',[1,1,1],'MarkerSize',4.5,'linewidth',1};
set_figsize(4,8*2,8)
subaxis(1,2,1,'S',.1,'ML',.0925,'MR',0.01,'MB',.135)
for ll = 1:numel(L)
    for ee = 1:Emax
        plot(minmax(B),[1,1]*real(eigs0{ll}(ee)),':','color',cols{4}); hold on
        plot(B,real(eigs0{ll}(ee) + delta{ll}(ee,:)),'-','color',Lcols{ll},'linewidth',1)
        plot(B,real(sorteigs(ee,:,ll)),markerstyle{:},'color',Lcols{ll})        
    end
end
patch(xlims([1,end,end,1]),e*[0,xlims([1,end]),0]*vf^2/(ef*ev2jo)*hbar_eV,cols{8})
hold off
ylim(ylimsRe); xlim(xlims)
xlabel('Magnetic field \itB\rm (T)','Fontsize',10)
ylabel('Re(\omega) (eV)','Fontsize',10)
set(gca,'YTick',0:.02:.1,'Layer','Top','LineWidth',.3,'XColor',cols{4},'YColor',cols{4})

%---- manual fixing of peculiar L=400 nm case
oldsorteigs = sorteigs;
sorteigs(2,[end-1,end],3) = real(sorteigs(2,[end-1,end],3)) + 1i*50; 
%---- fixing done 
subaxis(1,2,2)
for ll = 1:numel(L)
    for ee = 1:2
        plot(minmax(B),-[1,1]*imag(eigs0{ll}(ee)),':','color',cols{4}); hold on
        pApprox = plot(B,-imag(eigs0{ll}(ee) + delta{ll}(ee,:)),'-','color',Lcols{ll},'linewidth',1);
        pNumer = plot(B,-imag(sorteigs(ee,:,ll)),markerstyle{:},'color',Lcols{ll});
    end
end
hold off
ylim([4.1e-4,6e-4]); xlim(minmax(B)+[-1,1]*.25)
legend([pApprox,pNumer],{'\alpha(1-i\gamma/2\omega^0)\omega_c',...
    'Numerics'},'Location','NorthWest'); legend boxoff
set(gca,'YTick',(4.2:.4:5.8)*1e-4)
xlabel('Magnetic field \itB\rm (T)','Fontsize',10)
ylabel('-Im(\omega) (eV)','Fontsize',10,'LineWidth',.3,'Color',cols{4})


export_fig('Disk_DipoleSplittingUnderMagneticField','-pdf')