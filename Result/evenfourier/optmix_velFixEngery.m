function [u,v,hmix]=optmix_velFixEngery(that,Diff,Lx,Ly)%得到最优速度场，that：t的傅立叶变换
% compute the velocity for optimal mixing
% input: that the fft2 of scalar field
% n,m,kmag are the matrices related to fft2.
% output: u,v are velocity field in physical domain. 

%% Phi=Del^{-1}\theta in Fourier
M=size(that,1);
N=size(that,2);
Invthat=-that./Diff.kmag.^2;
Invthat(1,1)=0;

%% \Grad\phi
Invthatx=1i*Diff.m.*Invthat;%% Temp=phi k   Fourier space
Invthaty=1i*Diff.n.*Invthat;

hmix=sqrt(sum(sum(abs (Invthatx).^2+abs (Invthaty).^2 ,1),2))/N/M;
%% \theta\grad\phi in Real and Fourier (De-aliased for the
%% multiplication)
theta=real(ifft2(that));                  %% Temp=theta k，虚部只影响数值误差
u=fft2(theta.*real(ifft2(Invthatx)));
v=fft2(theta.*real(ifft2(Invthaty)));
    u=dealiasingf(u);
    v=dealiasingf(v);


%
%         unp2nm=sqrt(sum(sum(abs (u).^2+abs (v).^2,1),2))*2*pi/N^2;

%% Incomp. Proj. of Del^{-1} \theta\grad\phi in Fourier and in Real  投影算子
uhat=(ones(size(that))-Diff.m.*Diff.m./Diff.kmag.^2).*u-Diff.m.*Diff.n./Diff.kmag.^2.*v;
vhat=(ones(size(that))-Diff.n.*Diff.n./Diff.kmag.^2).*v-Diff.n.*Diff.m./Diff.kmag.^2.*u;
uhat(1,1)=0;
vhat(1,1)=0;
u=real (ifft2(uhat));
v=real (ifft2(vhat));
%the scale u^2+v^2 = Lx^2*Ly^2
rec=sqrt(f2int(uhat,Lx,Ly)+f2int(vhat,Lx,Ly));
u=u/rec;
v=v/rec;

end