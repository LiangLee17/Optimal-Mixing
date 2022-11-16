function y = f2int(fhat,Lx,Ly)
%input: the fft2 of the scalar field f
%rectanglar domain [0,Lx]*[0,Ly]

%output: the integral of f^2 over the rectanglar domain. 


if ~exist('Lx', 'var')
    Lx=2*pi;
end
if ~exist('Ly', 'var')
    Ly=2*pi;
end

%sum(sum(abs (fhat/N/M).^2,1),2)*pi^2 is the int f^2
% over [0,2*pi]*[0,2*pi]
% the factor Lx*Ly/(2*pi)^2
M=size(fhat,1);
N=size(fhat,2);
y=sum(sum(abs (fhat/N/M).^2,1),2)*Lx*Ly;
end

