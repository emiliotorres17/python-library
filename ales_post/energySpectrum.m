close all; 
clear all 
clc; 
format compact 
%#ok<*NOPTS>;

N = 128; 

name    = '../media/dynamic-smag.png';
Title   = 'Dynamic Smagorinsky $\nu=0.001324$'
lab     = '$\nu = 0.001324$'

% Loading data
u1  = load('../data/velocity-1.mat')
u2  = load('../data/velocity-2.mat')
u3  = load('../data/velocity-3.mat')
u1  = u1.u1(:,:,:,end);
u2  = u2.u2(:,:,:,end);
u3  = u3.u3(:,:,:,end);

% Input your velocity vector field data into the matrix "U"
U           = zeros([3,N,N,N]);
U(1,:,:,:)  = u1;
U(2,:,:,:)  = u2;
U(3,:,:,:)  = u3;

nx = [N N N];
nxc = [N/2+1 N/2+1 N/2+1];
k = 1:N/2+1;

%% Calculating Energy Spectrum
[Eu,Ev,Ew]=energySpectra(U,nx,nxc,N); 
Ek = Eu+Ev+Ew;
save('../spectra-N-128-DS-40-nu-0013.mat', 'Ek')
%% Plotting Guidelines
figure; hold on
plot(k(2:N/2),Ek(2:N/2)/Ek(2),'r', 'linewidth',2)
plot(k(2:N/2-12),10^1*k(2:N/2-12).^(-5/3),'k-.','linewidth',2)
ylim([10^-6.0, 1.5])
title(Title, 'interpreter', 'latex')
ylabel(strcat('$E(k)$'),'interpreter','latex','fontsize',17); 
xlabel(strcat('$k$'),'interpreter','latex','fontsize',17);
set(gca, 'XScale', 'log', 'YScale', 'log');
grid('on')
saveas(gca, name)

%% Function
function [Eu,Ev,Ew]=energySpectra(var,nx,nxc,N)
 
kmax=round(sqrt((1-nxc(1))^2+(1-nxc(2))^2+(1-nxc(3))^2))+1;
Eu=zeros(kmax,1);
Ev=zeros(kmax,1);
Ew=zeros(kmax,1);
 
%disp('Calculating Eu');
P=abs(fftshift(fftn(squeeze(var(1,:,:,:))))).^2;
for i=1:nx(1)
    for j=1:nx(2)
        for k=1:nx(3)
            km=round(sqrt((i-nxc(1))^2+(j-nxc(2))^2+(k-nxc(3))^2))+1;
            Eu(km)=Eu(km)+P(i,j,k);
        end
    end
end
 
%disp('Calculating Ev');
P   = abs(fftshift(fftn(squeeze(var(2,:,:,:))))).^2;
for i=1:nx(1)
    for j=1:nx(2)
        for k=1:nx(3)
            km=round(sqrt((i-nxc(1))^2+(j-nxc(2))^2+(k-nxc(3))^2))+1;
            Ev(km)=Ev(km)+P(i,j,k);
        end
    end
end
 
%disp('Calculating Ew');
P=abs(fftshift(fftn(squeeze(var(3,:,:,:))))).^2;
for i=1:nx(1)
    for j=1:nx(2)
        for k=1:nx(3)
            km=round(sqrt((i-nxc(1))^2+(j-nxc(2))^2+(k-nxc(3))^2))+1;
            Ew(km)=Ew(km)+P(i,j,k);
        end
    end
end
clear P;
 
end
