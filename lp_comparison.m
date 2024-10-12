%workspace of correlation code should be given as input 
% this code finds the Localization precision by merging the psf of
% fortunate molecules.
% 


stop=size(final_two,1);
wbox=13;
failed2=ones(stop,1);
for i=1:stop
    a=final_two(i,1);
    b=final_two(i,2);
    g1=grab_mol(:,:,a);
    g2=grab_mol(:,:,b);
    g_two=g1+g2;
    xc_box=((2*G_arr(a,17))+11)*0.5;
    yc_box=((2*G_arr(a,18))+11)*0.5;

        xguess=floor(G_arr(a,3))-xc_box;
        yguess=floor(G_arr(a,4))-yc_box;
        
        for i1=1:wbox
            for j1=1:wbox
                kz=(i1-1)*wbox+j1;
                xymerge(kz)=0;
                zmerge(kz)=g_two(i1,j1);
            end
        end 

        beta0=[xguess,yguess,50,psf_w0/q,min(g_two(:))]; 
[betafit,resid,J,COVB,mse] = nlinfit(xymerge,zmerge,@gaussian_merge3,beta0);
      yf=betafit(2)+yc_box;
        xf=betafit(1)+xc_box;
        a0=betafit(3);
        r0=abs(betafit(4));
        off=betafit(5);
        
       
        if(a0 < 0)
           a0=0;
            failed2(i,1)=0;
             xf=0;
            yf=0;
        end
        if(xf > G_arr(a,17)+11 || xf < G_arr(a,17) || yf > G_arr(a,18)+11 || yf < G_arr(a,18)) %this is also checking for failure but where it checked a0 before it checks to see if the
           
            failed2(i,1)=0;
            xf=0;
            yf=0;
            a0=0;
        end
        
        
        yf_new(i,1)=yf;
        xf_new(i,1)=xf;
        r0_new(i,1)=r0;
        a0_new(i,1)=a0;
        
       
end
TK=a0_new(:,1)==0 ;
a0_new(TK,:)=[]; 
xf_new(TK,:)=[];
yf_new(TK,:)=[];
r0_new(TK,:)=[];
npix_all=pi*(r0_new.^2);    % area of molecule in square pixels
N=npix_all.*a0_new; % number of photons for each molecule

%Calculating the parameters after UPDATE:|
            lp2_new=((r0_new*q).^2+(q^2)/12)*1./N+8*pi*((r0_new*q).^4)*(bkgn^2)/(q^2)*1./(N.*N);
            lp_new=1.3*sqrt(lp2_new); %loc prec in um (add 30%, Thompson et al, 2002) -- localization precision?
   G_two_new_7=[xf_new yf_new a0_new r0_new lp_new];
   lp_mean2=mean(lp_new)*1000;
%clearvars -except G_two_new TK


stop=size(final_three,1);
failed3=ones(stop,1);

for j=1:stop
    a=final_three(j,1);
    b=final_three(j,2);
    c=final_three(j,3);
    g1=grab_mol(:,:,a);
    g2=grab_mol(:,:,b);
    g3=grab_mol(:,:,c);
    g_three=g1+g2+g3;
    xc_box=((2*G_arr(a,17))+11)*0.5;
    yc_box=((2*G_arr(a,18))+11)*0.5;

        xguess=floor(G_arr(a,3))-xc_box;
        yguess=floor(G_arr(a,4))-yc_box;
        
        for i1=1:wbox
            for j1=1:wbox
                kz=(i1-1)*wbox+j1;
                xymerge(kz)=0;
                zmerge(kz)=g_three(i1,j1);
            end
        end 

        beta0=[xguess,yguess,50,psf_w0/q,min(g_three(:))]; 
[betafit,resid,J,COVB,mse] = nlinfit(xymerge,zmerge,@gaussian_merge3,beta0);
      yf=betafit(2)+yc_box;
        xf=betafit(1)+xc_box;
        a0=betafit(3);
        r0=abs(betafit(4));
        off=betafit(5);
        
       
        if(a0 < 0)
           a0=0;
            failed3(i,1)=0;
             xf=0;
            yf=0;
        end
        if(xf > G_arr(a,17)+11 || xf < G_arr(a,17) || yf > G_arr(a,18)+11 || yf < G_arr(a,18)) %this is also checking for failure but where it checked a0 before it checks to see if the
           
            failed3(i,1)=0;
            xf=0;
            yf=0;
            a0=0;
        end
        
        
        yf_new3(j,1)=yf;
        xf_new3(j,1)=xf;
        r0_new3(j,1)=r0;
        a0_new3(j,1)=a0;
        
        
        
end
TK3=a0_new3(:,1)==0 ;
a0_new3(TK3,:)=[]; 
xf_new3(TK3,:)=[];
yf_new3(TK3,:)=[];
r0_new3(TK3,:)=[];
npix_all3=pi*(r0_new3.^2);    % area of molecule in square pixels
N3=npix_all3.*a0_new3; % number of photons for each molecule

%Calculating the parameters after UPDATE:|
            lp2_new3=((r0_new3*q).^2+(q^2)/12)*1./N3+8*pi*((r0_new3*q).^4)*(bkgn^2)/(q^2)*1./(N3.*N3);
            lp_new3=1.3*sqrt(lp2_new3); %loc prec in um (add 30%, Thompson et al, 2002) -- localization precision?
   G_three_7=[xf_new3 yf_new3 a0_new3 r0_new3 lp_new3];
   lp_mean3=mean(lp_new3)*1000;


stop=size(final_four,1);
failed4=ones(stop,1);

for k=1:stop
    a=final_four(k,1);
    b=final_four(k,2);
    c=final_four(k,3);
    d=final_four(k,4);
    g1=grab_mol(:,:,a);
    g2=grab_mol(:,:,b);
    g3=grab_mol(:,:,c);
    g4=grab_mol(:,:,d);
    g_four=g1+g2+g3+g4;
    xc_box=((2*G_arr(a,17))+11)*0.5;
    yc_box=((2*G_arr(a,18))+11)*0.5;

        xguess=floor(G_arr(a,3))-xc_box;
        yguess=floor(G_arr(a,4))-yc_box;
        
        for i1=1:wbox
            for j1=1:wbox
                kz=(i1-1)*wbox+j1;
                xymerge(kz)=0;
                zmerge(kz)=g_four(i1,j1);
            end
        end 

        beta0=[xguess,yguess,50,psf_w0/q,min(g_four(:))]; 
[betafit,resid,J,COVB,mse] = nlinfit(xymerge,zmerge,@gaussian_merge3,beta0);
      yf=betafit(2)+yc_box;
        xf=betafit(1)+xc_box;
        a0=betafit(3);
        r0=abs(betafit(4));
        off=betafit(5);
        
       
        if(a0 < 0)
           a0=0;
            failed4(k,1)=0;
             xf=0;
            yf=0;
        end
        if(xf > G_arr(a,17)+11 || xf < G_arr(a,17) || yf > G_arr(a,18)+11 || yf < G_arr(a,18)) %this is also checking for failure but where it checked a0 before it checks to see if the
           
            failed4(k,1)=0;
            xf=0;
            yf=0;
            a0=0;
        end
        
        
        yf_new4(k,1)=yf;
        xf_new4(k,1)=xf;
        r0_new4(k,1)=r0;
        a0_new4(k,1)=a0;
        
        
        
end



TK4=a0_new4(:,1)==0 ;
a0_new4(TK4,:)=[]; 
xf_new4(TK4,:)=[];
yf_new4(TK4,:)=[];
r0_new4(TK4,:)=[];
npix_all4=pi*(r0_new4.^2);    % area of molecule in square pixels
N4=npix_all4.*a0_new4; % number of photons for each molecule

%Calculating the parameters after UPDATE:|
            lp2_new4=((r0_new4*q).^2+(q^2)/12)*1./N4+8*pi*((r0_new4*q).^4)*(bkgn^2)/(q^2)*1./(N4.*N4);
            lp_new4=1.3*sqrt(lp2_new4); %loc prec in um (add 30%, Thompson et al, 2002) -- localization precision?
   G_four_7=[xf_new4 yf_new4 a0_new4 r0_new4 lp_new4];
   lp_mean4=mean(lp_new4)*1000;
%clearvars -except G_two_new TK
% save('G_new_0.8');


G_arr1=[G_two_new_7;G_three_7;G_four_7];
save('G_arr_0.7.mat','G_arr1'); % specify the correlation factor while saving

