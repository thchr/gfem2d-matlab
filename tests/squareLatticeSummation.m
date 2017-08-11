clc; clear all; close all;
NL = 10:50:1500;
x = .2;
y = .1;
kx = 2*pi/3;
ky = 2*pi/3;

V = zeros(size(NL));
for cc = 1:numel(NL)
    N = NL(cc);
    for m = -N:N
        for n = -N:N
            V(cc) = V(cc) + exp(1i*kx*n)*exp(1i*ky*m)/sqrt((x+n)^2+(y+m)^2);
        end
    end
end

% V = zeros(size(NL));
% for cc = 1:numel(NL)
%     N = NL(cc);
%     for n = -N:N
%         for m = -N:N
%             V(cc) = V(cc) + 1/sqrt((x+n)^2+(y+m)^2);
%         end
%     end
% end


plot(NL,V)