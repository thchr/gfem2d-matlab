clear; close all; clc; cols=flatcolors;

Ns = 230; Ncirc = 41; adivd = 1.5*sqrt(3)*5/2;

% 
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
hole = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----
R = {[1,0],[cosd(60),sind(60)]};
expandshrink = 0;

inclusions=cell(1,6);
cntrs = [cos(pi/3*(1:6).'),sin(pi/3*(1:6).')]/3*(1+expandshrink);
for ii = 1:6
    inclusions{ii} = bsxfun(@plus, hole, cntrs(ii,:));
end

%% MESHING
fprintf('\n|----- MESHING UNIT CELL -----|\n')
[mesh,~,bnd] = geomPeriodicInclusion('triangular',inclusions,Ns,...
    @(x,y) 7/Ns*ones(size(x)),1);

%%
repetitions = -1:1;

try close(2); figure(2); catch; figure(2); end
hold on;
for rr1 = repetitions
    for rr2 = repetitions
        
        for ii = 1:numel(inclusions)
            plot(inclusions{ii}([1:end,1],1)+R{1}(1)*rr1+R{2}(1)*rr2,...
                inclusions{ii}([1:end,1],2)+R{1}(2)*rr1+R{2}(2)*rr2,...
                'Color',cols{4}*(rr1==0 & rr2==0) + cols{8}*(rr1~=0 | rr2~=0))
        end
        plot(cntrs(:,1),cntrs(:,2),'ok','MarkerFaceColor',cols{14});
        for ii = 1:6
            text(cntrs(ii,1),cntrs(ii,2),num2str(ii))
            
            plot(bnd([1:end,1],1)+R{1}(1)*rr1+R{2}(1)*rr2,...
                bnd([1:end,1],2)+R{1}(2)*rr1+R{2}(2)*rr2,'k-')
        end
    end
end

hold off;

axis equal
set(gca,'Fontsize',8)
xlim((minmax(repetitions)+[-1,1] + [-.075,.075])*abs(R{1}(1)+R{2}(1))/2)
ylim((minmax(repetitions)+[-1,1] + [-.075,.075])*abs(R{1}(2)+R{2}(2))/2)

