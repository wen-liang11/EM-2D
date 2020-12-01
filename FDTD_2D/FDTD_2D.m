clear, clc;

% choose 1st or 2nd order ABC to use
first_order=0;
second_order=0;
% choose source
gaussian=1;
impulse=0;

xsize=100;
ysize=100;
timetot=300;
xcenter=50;
ycenter=50;

epsilon=((1/(36*pi))*1e-9)*ones(xsize,ysize);     % permittivity
mu=(4*pi*1e-7)*ones(xsize,ysize);     % permeability
c=3e+8;

S=0.5;
delta=1e-6;     % dx = dy     
deltat=S*delta/c;

Ez=zeros(xsize,ysize);
Hy=zeros(xsize,ysize);
Hx=zeros(xsize,ysize);

sigma=4e-4*ones(xsize,ysize);       % electric conductivity
sigma_star=4e-4*ones(xsize,ysize);      % magnetic conductivity

A=((mu-0.5*deltat*sigma_star)./(mu+0.5*deltat*sigma_star));     % Coefficient in H matrix 
B=(deltat/delta)./(mu+0.5*deltat*sigma_star);
                                                    
C=((epsilon-0.5*deltat*sigma)./(epsilon+0.5*deltat*sigma));     % Coefficient in E matrix 
D=(deltat/delta)./(epsilon+0.5*deltat*sigma);  

c1=(c*deltat-delta)/(c*deltat+delta);     % Coefficient in ABC
c2=(4*delta^3-2*delta*(c*deltat)^2)/(2*delta^2*(delta+c*deltat));
c3=(delta*(c*deltat)^2)/(2*delta^2*(delta+c*deltat));

last_x_for=zeros(1,ysize);          % previous timesteps in ABC
last_x_minus1_for=zeros(1,ysize);
last_y_for=zeros(xsize,1);
last_y_minus1_for=zeros(xsize,1);
last_x_bac=zeros(1,ysize);
last_x_plus1_bac=zeros(1,ysize);
last_y_bac=zeros(xsize,1);
last_y_plus1_bac=zeros(xsize,1);

for n=1:1:timetot
     
    if gaussian==0 && n==1          % impulse input
        Ez(xcenter,ycenter)=1;
    end
    
    n1=2;
    n2=xsize-2;
    n11=2;
    n22=ysize-2;
    
    % Hy and Hx fields
    Hx(n1:n2-1,n11:n22-1)=A(n1:n2-1,n11:n22-1).*Hx(n1:n2-1,n11:n22-1)-B(n1:n2-1,n11:n22-1).*(Ez(n1:n2-1,n11+1:n22)-Ez(n1:n2-1,n11:n22-1));
    Hy(n1:n2-1,n11:n22-1)=A(n1:n2-1,n11:n22-1).*Hy(n1:n2-1,n11:n22-1)+B(n1:n2-1,n11:n22-1).*(Ez(n1+1:n2,n11:n22-1)-Ez(n1:n2-1,n11:n22-1));
    if first_order == 0 && second_order == 0
        Hx(:,2)=0;
        Hx(:,97)=0;
        Hy(2,:)=0;
        Hy(97,:)=0;
    end
 
    % Ez field
    Ez(n1+1:n2-1,n11+1:n22-1)=C(n1+1:n2-1,n11+1:n22-1).*Ez(n1+1:n2-1,n11+1:n22-1)+(Hy(n1+1:n2-1,n11+1:n22-1)-Hy(n1:n2-2,n11+1:n22-1)-Hx(n1+1:n2-1,n11+1:n22-1)+Hx(n1+1:n2-1,n11:n22-2)).*D(n1+1:n2-1,n11+1:n22-1);
    
    
    %1st order ABC
    if first_order == 1
        %forward boundary
        if n>=xsize-2-xcenter
            Ez(xsize-2,3:ysize-3)=c1*(Ez(xsize-3,3:ysize-3)-last_x_for(1,3:ysize-3))+last_x_minus1_for(1,3:ysize-3);
        end
        last_x_for(1,1:1:ysize)=Ez(xsize-2,1:1:ysize);
        last_x_minus1_for(1,1:1:ysize)=Ez(xsize-3,1:1:ysize);
        %backward boundary
        if n>=xcenter-3
            Ez(2,3:ysize-3)=c1*(Ez(3,3:ysize-3)-last_x_bac(1,3:ysize-3))+last_x_plus1_bac(1,3:ysize-3);
        end
        last_x_bac(1,1:1:ysize)=Ez(2,1:1:ysize);
        last_x_plus1_bac(1,1:1:ysize)=Ez(3,1:1:ysize);
        %upward boundary
        if n>=ysize-2-ycenter
            Ez(3:xsize-3,ysize-2)=c1*(Ez(3:xsize-3,ysize-3)-last_y_for(3:xsize-3,1))+last_y_minus1_for(3:xsize-3,1);
        end
        last_y_for(1:1:xsize,1)=Ez(1:1:xsize,ysize-2);
        last_y_minus1_for(1:1:xsize,1)=Ez(1:1:xsize,ysize-3);
        %downward boundary
        if n>=ycenter-3
            Ez(3:xsize-3,2)=c1*(Ez(3:xsize-3,3)-last_y_bac(3:xsize-3,1))+last_y_plus1_bac(3:xsize-3,1);
        end
        last_y_bac(1:1:xsize,1)=Ez(1:1:xsize,2);
        last_y_plus1_bac(1:1:xsize,1)=Ez(1:1:xsize,3);
    end
    
    
    % 2nd order ABC
    if second_order == 1
        %forward boundary
        if n>=xsize-2-xcenter
            Ez(xsize-2,3:1:ysize-3)=c1*(Ez(xsize-3,3:1:ysize-3)+last_last_x_for(1,3:1:ysize-3))-last_last_x_minus1_for(1,3:1:ysize-3)+c2*(last_x_for(1,3:1:ysize-3)+last_x_minus1_for(1,3:1:ysize-3))+c3*(last_x_minus1_for(1,2:1:ysize-4)+last_x_minus1_for(1,4:1:ysize-2)+last_x_for(1,2:1:ysize-4)+last_x_for(1,4:1:ysize-2));
        end

        %previous time steps updated at forward boundary
        last_last_x_for=last_x_for;
        last_last_x_minus1_for=last_x_minus1_for;
        last_x_for(1,1:1:ysize)=Ez(xsize-2,1:1:ysize);
        last_x_minus1_for(1,1:1:ysize)=Ez(xsize-3,1:1:ysize);

        %backward boundary
        if n>=xcenter-3
            Ez(2,3:1:ysize-3)=-last_last_x_bac(1,3:1:ysize-3)+c1*(Ez(3,3:1:ysize-3)+last_last_x_plus1_bac(1,3:1:ysize-3))+c2*(last_x_bac(1,3:1:ysize-3)+last_x_plus1_bac(1,3:1:ysize-3))+c3*(last_x_plus1_bac(1,2:1:ysize-4)+last_x_plus1_bac(1,4:1:ysize-2)+last_x_bac(1,2:1:ysize-4)+last_x_bac(1,4:1:ysize-2));
        end

        %previous time steps updated at backward boundary
        last_last_x_bac=last_x_bac;
        last_last_x_plus1_bac=last_x_plus1_bac;
        last_x_bac(1,1:1:ysize)=Ez(3,1:1:ysize);
        last_x_plus1_bac(1,1:1:ysize)=Ez(2,1:1:ysize);

        %upward boundary
        if n>=ysize-2-ycenter
            Ez(3:1:xsize-3,ysize-2)=c1*(Ez(3:1:xsize-3,ysize-3)+last_last_y_for(3:1:xsize-3,1))-last_last_y_minus1_for(3:1:xsize-3,1)+c2*(last_y_for(3:1:xsize-3,1)+last_y_minus1_for(3:1:xsize-3,1))+c3*(last_y_minus1_for(2:1:xsize-4,1)+last_y_minus1_for(4:1:xsize-2,1)+last_y_for(2:1:xsize-4,1)+last_y_for(4:1:xsize-2,1));
        end

        %previous time steps updated at upward boundary
        last_last_y_for=last_y_for;
        last_last_y_minus1_for=last_y_minus1_for;
        last_y_for(1:1:xsize,1)=Ez(1:1:xsize,ysize-2);
        last_y_minus1_for(1:1:xsize,1)=Ez(1:1:xsize,ysize-3);

        %downward boundary
        if n>=ycenter-3
            Ez(3:1:xsize-3,2)=-last_last_y_bac(3:1:xsize-3,1)+c1*(Ez(3:1:xsize-3,3)+last_last_y_plus1_bac(3:1:xsize-3,1))+c2*(last_y_bac(3:1:xsize-3,1)+last_y_plus1_bac(3:1:xsize-3,1))+c3*(last_y_plus1_bac(2:1:xsize-4,1)+last_y_plus1_bac(4:1:xsize-2,1)+last_y_bac(2:1:xsize-4,1)+last_y_bac(4:1:xsize-2,1));
        end

        %previous time steps updated at downward boundary
        last_last_y_bac=last_y_bac;
        last_last_y_plus1_bac=last_y_plus1_bac;
        last_y_bac(1:1:xsize,1)=Ez(1:1:xsize,3);
        last_y_plus1_bac(1:1:xsize,1)=Ez(1:1:xsize,2);

        %corner update
        Ez(2,2)=last_last_x_bac(3);
        Ez(2,ysize-2)=last_last_x_bac(ysize-3);
        Ez(xsize-2,2)=last_last_x_minus1_for(3);
        Ez(xsize-2,ysize-2)=last_last_x_minus1_for(ysize-3);
    end
    
    % gaussian
    if gaussian==1
        if n<=42
            Ez(xcenter,ycenter)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/32;
        else
            Ez(xcenter,ycenter)=0;
        end
    end
    
        % impulse
    if impulse == 1 && n == 1
        Ez(xcenter,ycenter)=0;
    end
    
    %plot of Ez
    s1= surf(delta*(1:1:xsize)*1e+6,(1e+6*delta*(1:1:ysize))',Ez);
    s1.EdgeColor='none';
    title(['\fontsize{12}Ez vs X for 2D FDTD & 2nd ABC at time = ',num2str(round(n*deltat*1e+15)),' fs']); 
    xlabel('x (in um)','FontSize',12);
    ylabel('y (in um)','FontSize',12);
    set(gca,'FontSize',14,'color','w');
    axis([0 100 0 100 -0.5 0.5]);
    F(n)=getframe(gcf);
end
% writerObj=VideoWriter('test1.avi');
% open(writerObj);
% writeVideo(writerObj,F);
% close(writerObj);