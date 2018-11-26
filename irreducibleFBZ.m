function [kvec,kplot,kmark,FBZ,irrFBZ] = irreducibleFBZ(type,Nmin,plotchoice)
%CALL: [kvec,kplot,kmark,FBZ,irrFBZ] = irreducibleFBZ(type,Nmin,plotchoice)
%DESCRIPTION: Provides a list of k-points in the irreducible Brillouin zone 
%             of the specified type of lattice ('type') at Nmin kpoints.
%INPUT:  type | indicate type of lattice (assumed in unit-form). Options:
%               'square': Square lattice of unit-side
%                         Matches R{1}=[1,0] & R{2}=[0,1]
%               'triangular': Triangular lattice of unit side with
%                             hexagonal Brillouin zone (Gamma, M, K)
%                             Matches R{1}=[1,0] & R{2}=[cosd(60),sind(60)]
%        Nmin | Minimum number of points along the irreducible Brillouin
%               zone. Output may exceed this number of points by a few, in
%               certain nefarious cases.
%        plotchoice | Set to '1' or string 'plotfbz' to show the (irr)FBZ.
%OUTPUT: kvec | A list of (approx) Nmin xy-coordinates [size: Nmin x 2]
%               along the irreducible Brillouin zone. Start and end points
%               will coincide to ensure a full loop; otherwise no
%               duplicates.
%        kplot | A column vector of lengths [size: Nmin x 1] corresponding
%                to the length traversed around the irreducible Brillouin
%                zone (except for final point, which is padded to avoid
%                duplicates in kplot). 
%        kmark | A structure with fields:
%                symbol | A cell array of the symbols corresponding to
%                         high-symmetry points in the irreducible FBZ path.
%                n | A vector array indicating the position of high-
%                    symmetry points (in kvec) in its first column and 
%                    their associated symbol (in kmark.symbol) in its
%                    second column.
%        FBZ | High symmetry points in the Brillouin zone.
%        irrFBZ | High symmetry points in the irreducible Brillouin zone.
%AUTHOR:                                THOMAS CHRISTENSEN (21 March, 2016)


if ~exist('type','var') || isempty(type) %Default value for type is a square
    type = 'square';
end

if ~exist('Nmin','var') || isempty(Nmin) %Default (minimum) number of k-points
    Nmin = 40;
end

%% TYPE OF LATTICE 
switch type
    case 'square'
        FBZ = [-1,-1; 1,-1; 1,1; -1,1]*pi; %FBZ of unit-length square lattice
        irrFBZ = [0,0; 1,0; 1,1]*pi; %Irreducible FBZ high symmetry points
        kmark.symbol = {'\Gamma','X','M'};
    case 'triangular'
        FBZ = 2*pi*[-2/3, 0; -1/3, -1/sqrt(3); 1/3, -1/sqrt(3); 2/3, 0; 1/3, 1/sqrt(3); -1/3, 1/sqrt(3)];
        irrFBZ = 2*pi*[0,0; 1/2,-1/(2*sqrt(3)); 2/3,0];
        kmark.symbol = {'\Gamma','M','K'};
end
irrFBZcyc = circshift(irrFBZ,-1);
irrFBZsteplen = sqrt( (irrFBZ(:,1) - irrFBZcyc(:,1) ).^2 + (irrFBZ(:,2) - irrFBZcyc(:,2) ).^2 );
irrFBZtotlen = sum(sqrt( (irrFBZ(:,1) - irrFBZcyc(:,1) ).^2 + (irrFBZ(:,2) - irrFBZcyc(:,2) ).^2 ));

%% CONSTRUCTION OF k-VECTOR AND ITS RELATED QUANTITIES
Ninc = Nmin;
for ee=1:2 %Ensure that we distribute at least N points on the boundary
    Ncur = round(Ninc*irrFBZsteplen/irrFBZtotlen);
    if sum(Ncur) >= Nmin
        break
    else
        Ninc = Ninc+1;
    end
end
%Add extra points where needed
Ncur(2:end) = Ncur(2:end)+1;

%N = N; %Readjust number of points to account for removal of "double" points at vertices in irrFBZ

kvec = []; %Preallocate
kmark.n = [1, 1];
for kcur = 1:size(irrFBZ,1)
    if kcur == size(irrFBZ,1)
        knext = 1;
    else
        knext = kcur+1;
    end
    %Ncur = round(N*norm(irrFBZ(kcur,:)-irrFBZ(knext,:),2)/irrFBZtotlen);
    kvecadd = [linspace(irrFBZ(kcur,1),irrFBZ(knext,1),Ncur(kcur)).', ...
               linspace(irrFBZ(kcur,2),irrFBZ(knext,2),Ncur(kcur)).'];
    if kcur == 1
        kvec = [kvec; kvecadd(1:end,:)];
    else %Avoid double points at vertex points (apart from start and end)
        kvec = [kvec; kvecadd(2:end,:)];
    end
    kmark.n = [kmark.n; [size(kvec,1),knext]]; %The indices n show where in kvec the high-symmetry points are (first column, indexed into kplot) and what symbol they are associated with (second column, indexed into kmark.symbol).
end

%Construct a vector, whose value equals the length in the loop around the
%irreducible Brillouin zone
kplot = cumsum(sqrt((kvec(:,1)-circshift(kvec(:,1),-1)).^2+(kvec(:,2)-circshift(kvec(:,2),-1)).^2));
kplot(end) = 2*kplot(end-1) - kplot(end-2); %The last point is equal to the starting point, so their interdistance is zero: undesirable for plotting; fixed manually here
kplot = kplot - kplot(1); %Readjust so that kplot(1) = 0 and k(end) = irrFBZtotlen
%% PLOTTING
if exist('plotchoice','var') && (all(plotchoice == 1) || strcmpi(plotchoice,'plotfbz'))
    cols=flatcolors;
    
    set_figsize([],12,12);
    hold on
    patch([FBZ(:,1);FBZ(1,1)],[FBZ(:,2);FBZ(1,2)],cols{14}*.15+.85,'LineStyle',':','EdgeColor',cols{2}*.85+.15)
    patch(irrFBZ(:,1),irrFBZ(:,2),cols{14}*.35+.5,'LineStyle','-','EdgeColor',cols{14})
    plot(irrFBZ(:,1),irrFBZ(:,2),'o','MarkerFaceColor',cols{14}*.2+.8,'MarkerEdgeColor',cols{14},'MarkerSize',8,'LineWidth',1.25)
    plot(kvec(:,1),kvec(:,2),'.','color',cols{14},'MarkerSize',5);
    hold off
    
    axis equal
    xlim(minmax(FBZ(:,1))+max(abs(FBZ(:,1)))*.1*[-1,1])
    ylim(minmax(FBZ(:,2))+max(abs(FBZ(:,2)))*.1*[-1,1])
    box on
    set(gca,'Fontsize',8,'LineWidth',.2)
    
    drawnow
end