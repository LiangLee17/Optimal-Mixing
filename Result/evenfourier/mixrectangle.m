%mix in rectangluar domain.
clear
close all

%%
%spatial parameters
M =256+1;
N =256+1;
ME=2*M-2;
NE=2*N-2;
Lx=2;
Ly=2;
x0=-1;
y0=-1;
x=linspace(x0,x0+Lx,M);
y=linspace(y0,y0+Ly,N);
[yy,xx]=meshgrid(y,x);

[Diff.n,Diff.m]=meshgrid(-ME/2:ME/2-1,-NE/2:NE/2-1);
Diff.m=pi/Lx*fftshift(fftshift(Diff.m,1),2);
Diff.n=pi/Ly*fftshift(fftshift(Diff.n,1),2);
Diff.kmag=sqrt(Diff.m.^2+Diff.n.^2);

%%
%initial function
testcase=4;
switch testcase
    case 1
        theta=cos(2*pi*xx).*cos(pi*yy)+0.5*cos(0*pi*xx).*cos(2*pi*yy);
        energy=0.15*Lx*Ly;
    case 2
        theta=(cos(2*pi*xx).*cos(1*pi*yy)+0.5*cos(0*xx).*cos(3*pi*yy))+xx+yy;
        energy=0.3*Lx*Ly;
    case 3
        theta=(cos(2*xx*2*pi/Lx).*cos(yy*2*pi/Ly)+0.5*cos(0*xx*2*pi/Lx).*cos(3*yy*2*pi/Ly));
        energy=0.15*Lx*Ly;
    case 4
        theta=sin(2*pi*xx).*sin(pi*yy)+0.5*sin(0*pi*xx).*sin(2*pi*yy);
        energy=0.15*Lx*Ly;
end
% thetaE = [theta; flipud(theta(2:end-1,:))];
% thetaE = [thetaE fliplr(thetaE(:,2:end-1))];
thetaE = EvenExtension(theta,'e');
that=fft2(thetaE);
that(1,1)=0;%shfit theta, such that it is mean zeron function.
%%
%time parameters
dt=0.02;
tt=0:dt:1;


divu=zeros(1,length(tt));
nl2=zeros(1,length(tt));
nmix=zeros(1,length(tt));
nhm1=zeros(1,length(tt));



%%
%evolution
for indt=1:length(tt)
    tic
    nl2(indt)=sqrt(f2int(that,Lx,Ly));%average of l2 norm
    if nl2(indt)>1.05*nl2(1) %jump out of this loop if explosion
        %the system explodes when reach tt(i)
        display(['The system explodes at ' num2str(tt(indt)) ])
        indt=indt-1;%don't save the wrong result
        break;
    end
    
    
    thatx=1i*Diff.m.*that;% theta_x
    thaty=1i*Diff.n.*that;
    %     nl2(i)^2./sqrt(f2int(thatx,Lx,Ly)+f2int(thaty,Lx,Ly)+nl2(i)^2);
    %    is the standard definition of h^{-1} norm, but theta^2, thetax^2 have
    % different scale. we'd like something invariant under the scaling of
    % domain's length.
    nhm1(indt)=1./sqrt(f2int(thatx,Lx,Ly)+f2int(thaty,Lx,Ly));%h^{-1} norm
    
    [u,v,nmix(indt)]=optmix_velFixEngery(that,Diff,Lx,Ly);
    u=u*energy;%scale velocity
    v=v*energy;
    u([1,end/2+1],:)=0;
    v(:,[1,end/2+1])=0;
    
    
    data.theta{indt}=theta;
    data.u{indt}=u(1:M,1:N);
    data.v{indt}=v(1:M,1:N);
    
    
    [~,Xs]=ode45(@(t,y0)velocity(t,y0,u,v,Diff),tt(indt):dt/2:tt(indt)+dt,that(:));
    that=reshape(Xs(end,:),ME,NE);
    that=dealiasingf(that);
    that(1,1)=0;
    theta=EvenExtension(real(ifft2(that)),'c');
    toc
end
fname=strcat('DataFixEnergy','Case',num2str(testcase),'_Ngrid',num2str(N), ...
    '_T',num2str(tt(indt)),'_dt',num2str(dt),'.mat');
save(fname,'xx','yy','data');

%%
% %plot the norms
figure (3)
semilogy(tt(1:indt),nl2(1:indt)/nl2(1),'-g')
hold on
semilogy(tt(1:indt),nmix(1:indt)/nmix(1),'-.r')

semilogy(tt(1:indt),nhm1(1:indt)/nhm1(1),'--b')
hold off
axis([tt(1),tt(indt),0,1.05])
legend('relative variance','relative mix norm','relative hm1 norm')


%
% % title('Stokes Flow','fontsize',14);
% % legend('\beta=0 (3D potential flow)','\beta=0.5','fontsize',N2)
% % xlabel('$\frac{\lambda}{l}$','interpreter','latex','fontsize',N2);
fname=strcat('NormsFixEnergy','Case',num2str(testcase),'_Ngrid',num2str(N), ...
    '_T',num2str(tt(indt)),'_dt',num2str(dt),'.jpg');
saveas(3,fname)
%%

MovScalar=moviein(length(tt));

figure(1)
for indplot=1:indt
    pcolor(xx,yy,data.theta{indplot})% the horizontal is x.
    shading interp
    colormap gray;
    % imshow(theta)
    MovScalar(:,indplot)=getframe(1);
end



fname=strcat('ScalarFixEnergy','Case',num2str(testcase),'_Ngrid',num2str(N), ...
    '_T',num2str(tt(indt)),'_dt',num2str(dt),'.gif');
%  movie(M,1,1)
dt1=5*dt;
for indgif=1:length(MovScalar)
    [image,map] = frame2im(MovScalar(indgif));
    [im,map2]=rgb2ind(image,128);
    if indgif==1
        imwrite(im,map2,fname,'GIF','WriteMode','overwrite','DelayTime',dt1,'LoopCount',inf);
    else
        imwrite(im,map2,fname,'WriteMode','append','DelayTime',dt1);
    end
end
%%
MovVel=moviein(length(tt));
figure(2)
for indplot=1:indt
    quiver(xx(1:floor(M/50):end,1:floor(N/50):end), ...
        yy(1:floor(M/50):end,1:floor(N/50):end), ...
        data.u{indplot}(1:floor(M/50):end,1:floor(N/50):end), ...
        data.v{indplot}(1:floor(M/50):end,1:floor(N/50):end),1.75)
    axis tight
    MovVel(:,indplot)=getframe(2);
end

fname=strcat('VelFixEnergy','Case',num2str(testcase),'_Ngrid',num2str(N), ...
    '_T',num2str(tt(indt)),'_dt',num2str(dt),'.gif');

dt1=5*dt;
for indgif=1:length(MovVel)
    [image,map] = frame2im(MovVel(indgif));
    [im,map2]=rgb2ind(image,128);
    if indgif==1
        imwrite(im,map2,fname,'GIF','WriteMode','overwrite','DelayTime',dt1,'LoopCount',inf);
    else
        imwrite(im,map2,fname,'WriteMode','append','DelayTime',dt1);
    end
end