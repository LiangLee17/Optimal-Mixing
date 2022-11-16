function  dy = velocity(t,that,u,v,Diff)
M=size(u,1);
N=size(u,2);
that=reshape(that,M,N);
thatx=1i*Diff.m.*that;%% Temp=phi k   Fourier space
thaty=1i*Diff.n.*that;


dy=-fft2(ifft2(thatx).*u+ifft2(thaty).*v);


% kappa=1e-6;
% thatxx=sqrt(-1)*m.*thatx;%% Temp=phi k   Fourier space
% thatyy=sqrt(-1)*n.*thaty;
% dy=dy+kappa*(thatyy+thatxx);


dy=dealiasingf(dy);
dy=dy(:);
end

