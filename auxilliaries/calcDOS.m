function [dos,ene_eV] = calcDOS(eig_eV,ene_eV,A)
%CALL:           [dos,ene_eV] = calcDOS(eig_eV,ene_eV,A)
%DESCRIPTION: Calculates the density of states 'dos' for eigenvalues 
%(energies) supplied in 'eig_eV' at energies supplied in 'ene_eV'. 
%INPUT: The following dimensions are expected: 
%	eig_eV | rows indicate eigenvalues in eV for fixed parameter space,
%            while columns may indicate distinct parameter spaces, i.e.
%            size is [N x M] (N = #number of eigenvals and M = param vals)
%   ene_eV | column vector of energies [J x 1] in eV
%   A      | If a unit-cell volume/area/length 'A' is supplied (as scalar),
%             the output of 'dos' has the correct units of 1/[eV m^dimA]
%OUTPUT: dos | Density of states in [J x M] array, i.e. each column gives
%              the DOS at the mth (m = 1,..,M) parameter space at [with
%              energies ene_eV(j = 1,..,J) running down along rows]
%AUTHOR: Thomas Christensen                                  11 May, 2016


%----- DEFAULTS AND CHECK -----
if all(imag(eig_eV(:)) == 0)       %If the supplied eigenvalues 
    eig_eV = eig_eV - 1i*(6e-3)/2; %are purely real, we artificially 
end                                %add some broadening.

if ~exist('ene_eV','var') || isempty(ene_eV) %In case ene_eV is not specified
    maxloss_eV = max(abs(imag(eig_eV(:))));
    ene_eV = linspace(min(real(eig_eV(:))) - maxloss_eV,max(real(eig_eV(:))) + maxloss_eV,1000).';
end
if isrow(eig_eV); eig_eV = eig_eV.'; end %If eig_eV was erroneously supplied as row vector, we flip it

if nargin < 3 %In case A is not specified we simply specify it as unit
    A = 1;    %[the units of dos are then not meaningful]
end



%----- CALCULATE DOS -----
dos = zeros(numel(ene_eV),size(eig_eV,2)); %Preallocate space
for m = 1:size(eig_eV,2) %Calculate by looping (could be sped up by array'ing, but nvm)
    for n = 1:size(eig_eV,1)
        if ~isnan(eig_eV(n,m)); %If any of the eigenenergies are NaN we interpret this as an "missing" mode, and do not include a contribution from it
            dos(:,m) = dos(:,m) + imag(1./(-eig_eV(n,m)+ene_eV));
        end
    end
end
dos = -2/(pi*A)*dos; %Units of eV^(-1) m^(-dimA) with dimA = dimension of A [Divide by 10^(9*dimA) to get to eV^(-1) nm^(-dimA)]