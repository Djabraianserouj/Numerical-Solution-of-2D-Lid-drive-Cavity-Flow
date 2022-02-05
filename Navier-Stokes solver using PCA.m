%%% Navier-Stokes Solver using the Predictor-Corrector Approach (PCA) %%%
clear all 
close all
clc

Nx=11;               %nb of nodes in x-direction               
Ny=11;               %nb od nodes in y-direction
L=1;                 %Dimension of domain
T=0.1;               %Total time
dt=0.0001;           %Time step size  
nt=T/dt;             %Total number of timesteps 
Re=0.1;              %Reynold's number
iterations=50000;    %Maximum number of iterations for Pressure correction                   
dx=L/(Nx-1);
dy=L/(Ny-1);
[X,Y]=meshgrid(0:dx:L,0:dy:L);
u=zeros(Nx+1,Ny+2);
v=zeros(Nx+2,Ny+1);
p=zeros(Nx+2,Ny+2);
u_star=zeros(Nx+1,Ny+2);
v_star=zeros(Nx+2,Ny+1);
u_final=zeros(Nx+1,Ny+1);
v_final=zeros(Nx+1,Ny+1);
figure 1;
for k=1:nt           %Iterating over time 
  k
    if (rem(k,5)==0)
      %contourf(rot90(rot90(rot90(p))));    %Pressure Contour
      %contourf(flipud(rot90(uv_final)));   %Velocity Contour
      %colorbar;
      quiver(flipud(rot90(u_final)),flipud(rot90(v_final)),5);    %Streamlines
      drawnow;
    endif
##      if (k==2)
##        h1 = subplot( 321);
##        contourf(rot90(rot90(rot90(p))));
##        colorbar;
##        title (h1,'Pressure contour at t=0s');
##        xlabel(h1,'x'); 
##        ylabel(h1,'y');        
##      endif
##      if (k==200)
##        h2 = subplot( 322);
##        contourf(rot90(rot90(rot90(p))));
##        colorbar;
##        title (h2,'Pressure contour at t=0.02s');
##        xlabel(h2,'x'); 
##        ylabel(h2,'y');        
##      endif
##      if (k==400)
##        h3 = subplot( 323);
##        contourf(rot90(rot90(rot90(p))));
##        colorbar;
##        title (h3,'Pressure contour at t=0.04s');
##        xlabel(h3,'x'); 
##        ylabel(h3,'y');        
##      endif
##      if (k==600)
##        h4 = subplot( 324);
##        contourf(rot90(rot90(rot90(p))));
##        colorbar;
##        title (h4,'Pressure contour at t=0.06s');
##        xlabel(h4,'x'); 
##        ylabel(h4,'y');        
##      endif
##      if (k==800)
##        h5 = subplot( 325);
##        contourf(rot90(rot90(rot90(p))));
##        colorbar;
##        title (h5,'Pressure contour at t=0.08s');
##        xlabel(h5,'x'); 
##        ylabel(h5,'y');        
##      endif
##      if (k==1000)
##        h6 = subplot( 326);
##        contourf(rot90(rot90(rot90(p))));
##        colorbar;
##        title (h6,'Pressure contour at t=0.1s');
##        xlabel(h6,'x'); 
##        ylabel(h6,'y');        
##      endif
    
    %%% BC's for momentum equations %%%
    u(:,1)=-u(:,2);           %Bottom   
    u(:,Ny+2)=2-u(:,Ny+1);    %Lid
    v(1,:)=-v(2,:);           %Right
    v(Nx+2,:)=-v(Nx+1,:);     %Left
    %%% X-momentum equation %%%
    for i=2:Nx
      for j=2:Ny+1
        diffusion = (1/Re)*((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+(u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2);
        con_x = (-0.25)*((u(i+1,j)+u(i,j)).^2-(u(i,j)+u(i-1,j)).^2)/dx ;
        con_y = (-0.25)*((u(i,j+1)+u(i,j)).*(v(i+1,j)+v(i,j))-(u(i,j)+u(i,j-1)).*(v(i+1,j-1)+v(i,j-1)))/dy;
        u_star(i,j)=u(i,j)+dt*(diffusion + con_x + con_y);                      
      endfor
    endfor  
    
    %%% Y-momentum equation %%%
    for i=2:Nx+1
      for j=2:Ny
        diffusion = (1/Re)*((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+(v(i,j+1)+-2*v(i,j)+v(i,j-1))/dy^2);
        con_x = (-0.25)*((u(i,j+1)+u(i,j)).*(v(i+1,j)+v(i,j))-(u(i-1,j+1)+u(i-1,j)).*(v(i,j)+v(i-1,j)))/dx;
        con_y = (-0.25)*((v(i,j+1)+v(i,j)).^2-(v(i,j)+v(i,j-1)).^2)/dy;
        v_star(i,j)=v(i,j)+dt*(diffusion + con_x + con_y);             
      endfor
    endfor    
    pt=p;         %Temporary Pressure Matrix to minimize the error
    for it=1:iterations	               %Iterate until a desired error is reached
      %it
      %%% BC's for pressure %%%
      pt(1,:)=pt(2,:);
      pt(Nx+2,:)=pt(Nx+1,:);
      pt(:,1)=pt(:,2);
      pt(:,Ny+2)=pt(:,Ny+1);
      
      %%% Poisson's equation %%%
      i=2:Nx+1; j=2:Ny+1; 
          pt(i,j) = 0.25*(pt(i+1,j) + pt(i-1,j) + pt(i,j+1) + pt(i,j-1) - (dx/dt)*(u_star(i,j)-u_star(i-1,j)+v_star(i,j)-v_star(i,j-1)));      
      err=max(max(pt-p));
      p=pt;
      if err<10^-5
          break;   
      endif
    endfor
    %%% Velocities u and v at timestep (n+1) %%%
    u(2:Nx,2:Ny+1) = u_star(2:Nx,2:Ny+1)-(dt/dx)*(p(3:Nx+1,2:Ny+1)-p(2:Nx,2:Ny+1));
    v(2:Nx+1,2:Ny) = v_star(2:Nx+1,2:Ny)-(dt/dy)*(p(2:Nx+1,3:Ny+1)-p(2:Nx+1,2:Ny));
    
    %%% Velocity correction to our actual grid points %%%
    u_final(1:Nx+1,1:Ny+1)=0.5*(u(1:Nx+1,2:Ny+2)+u(1:Nx+1,1:Ny+1));
    v_final(1:Nx+1,1:Ny+1)=0.5*(v(2:Nx+2,1:Ny+1)+v(1:Nx+1,1:Ny+1)); 
    uv_final=sqrt(u_final.^2+v_final.^2);
endfor
