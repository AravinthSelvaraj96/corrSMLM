% This program detects single molecule spots after background subtraction
% and fits 2D gaussian then find centroids of the molecule, radius of the spot and 
% Number of photons detected on a spot. 



clear all
clc 
warning off
global xpix ypix wbox; %sets global variables that can be referenced in other m files

% provide info about where the raw file is and where the result should be saved
data_dir = '/Volumes/Aravinth/d2at_folder/'; %raw file folder address
an_dir =  '/Volumes/Aravinth/d2at_folder/result/'; %Output file will be put here
base_name='d2at'; %name of the tif file
imagefile = strcat(data_dir,base_name,'.tif'); 

rbox=6;

q=0.06; % effective pixel size of the camera interms of um
wvlnth=580/1000; %convert wavelength from nm to um -- where um is micrometers  
NA=1.45; %NA of the objective lens
psf_scale=1; %The scale of the point spread function

n_start =1; %The first frame 
n_end   =10000;%The last frame 

clims=[0 200]; %2D array indicating scaling of the image display


pix_to_pho = 18; %pixel to phton ratio of EMCCD camera used
min_thresh = 13; %minimum threshold for a "bright spot"    
box_overlap_factor = 1.5; 

w_mask = round(rbox*box_overlap_factor); 
wbox=2*rbox+1; 
[xpix0,ypix0] = meshgrid(-2*rbox:2*rbox,-2*rbox:2*rbox); 
[xpix,ypix] = meshgrid(-rbox:rbox,-rbox:rbox);

psf_w0 = psf_scale*0.61*wvlnth/NA; 
psf_std=psf_w0/2; 
psf_w02=(psf_w0/q)*(psf_w0/q); 

%rolling ball parameters 
rball=6; 
se = strel('ball',rball,rball,0);  this is a strel data type
FWHM=1; 
rk=(FWHM)/sqrt(2*log(2)); 
kw=20; 
[Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2); 
kd=sqrt(Xgs.*Xgs+Ygs.*Ygs);
gs=exp(-2*kd.*kd/(rk*rk)); 
gs=gs/sum(sum(gs)); 

%initialize arrays
n_init=500000; 
xcm_all=zeros(n_init,1); 
xf_all=zeros(n_init,1); %x coordinate of the molecules
yf_all=zeros(n_init,1); %y coordinate of the molecules
a0_all=zeros(n_init,1); % Amplitude of the fitted gaussian on molecules
r0_all=zeros(n_init,1); % radius of the spot
off_all=zeros(n_init,1); 
frame_num_all=zeros(n_init,1); % molecule belong to which frame
xf_err_all=zeros(n_init,1); 
yf_err_all=zeros(n_init,1); 
a0_err_all=zeros(n_init,1); 
r0_err_all=zeros(n_init,1); 
off_err_all=zeros(n_init,1); 
grab_sum_all=zeros(n_init,1); %sum of gray levels of the spot.
box_matrix_xtop=zeros(n_init,1);
box_matrix_ytop=zeros(n_init,1);%takes top right corner coordinates
total_molecules=0; %there are no molecules!
n_fail_a0=0; % maybe this is the number of failed curve fits.
n_fail_outbox=0; %Related to failure but not sure how.

count=0;
grab_mol=zeros((2*rbox+1),(2*rbox+1),50000); % to get the each spot as a matrix


% finding background noise

if(exist('bkgn','var')==0) 
    bkgn=bg_noise_calc01(imagefile,pix_to_pho); 
        
    answer = questdlg(sprintf(['Noise: ',num2str(bkgn,'%2.2f'),', Thresh: ',num2str(min_thresh),'. Yes to continue, No to restart.']));
       
    if ~strcmp(answer,'Yes') %if it IS okay
        return; %Return the calculated background to matlab
    end 
end 


tic %begins time for performance analysis.

wb_norm = n_end-n_start; 
if wb_norm == 0 
    wb_norm = 1; 
end 



for fileloop=n_start:n_end
    compl = (fileloop-n_start)/wb_norm; 
    
    drawnow; %Update the waitbar on the screen.
    
  
    i1=double(imread(imagefile,fileloop))/pix_to_pho; 
    
    [yw,xw]=size(i1); 
    
    y_cord=1:yw;
    x_cord=1:xw;

        
    i1_gs = uint16(conv2(i1,gs,'same')); 
    bkg = double(imopen(i1_gs,se)); 
    
    iprod=i1-bkg; %subtracts the background from the original image 
    iprod=iprod.*(iprod>0); %set negative values to 0 -- so if you end up subtracting more light than was there to begin with, just turn it black.
    

    
    high_pixel_mask = zeros(yw,xw); %initializes the high pixel mask matrix to zeroes
    n_boxes=0; 
    boxes_xy=zeros(10000,3); 
   
    
    % finding the spots based on threshold
    [iy,ix] = find(iprod>= min_thresh); %find all pixels in iprod above the threshold (returns the indecies of such pixels.)
   
    Aixy=[iy,ix];

    high_pix_inframe = size(iy,1); 
    high_pix_index = 1;
    while high_pix_index <= high_pix_inframe 
        pix_val = 0; 
        
        if (iprod(iy(high_pix_index),ix(high_pix_index))>pix_val && high_pixel_mask(iy(high_pix_index),ix(high_pix_index))==0  && iy(high_pix_index)<yw-rbox-1 && iy(high_pix_index)>rbox+1  && ix(high_pix_index)<xw-rbox-1 && ix(high_pix_index)>rbox+1)
            pix_val = iprod(iy(high_pix_index),ix(high_pix_index));
            high_pixel_y = iy(high_pix_index);  
            high_pixel_x = ix(high_pix_index);
           
            
            high_pix_index = high_pix_index + 1;
        else
            high_pix_index = high_pix_index + 1;
            continue
        end
        if pix_val < min_thresh 
            break 
        end 


        x0_box=high_pixel_x-rbox;
        y0_box=high_pixel_y-rbox;
        x1_box=high_pixel_x+rbox;
        y1_box=high_pixel_y+rbox;

        x0_mask=high_pixel_x-w_mask; 
        if x0_mask < 1 
            x0_mask = 1; 
        end 
        x1_mask=high_pixel_x+w_mask; 
        if x1_mask > xw
            x1_mask = xw;
        end
        y0_mask=high_pixel_y-w_mask; 
        if y0_mask < 1
            y0_mask = 1;
        end
        y1_mask=high_pixel_y+w_mask;
        if y1_mask > yw
            y1_mask = yw;
        end

        high_pixel_mask(y0_mask:y1_mask,x0_mask:x1_mask)=1; 

        grab=iprod(y0_box:y1_box,x0_box:x1_box); 
        grab_sum=sum(sum(grab)); 
      
        xm_sum=0;
        ym_sum=0;
        m_sum=0;

        for i=x0_box:x1_box 
            for j=y0_box:y1_box 
                xind=floor(i);
                yind=floor(j);
                intens=iprod(yind,xind); 
                xm_sum=xm_sum+xind*intens; 
                ym_sum=ym_sum+yind*intens;
                m_sum=m_sum+intens;
            end
        end

        x_cm=xm_sum/m_sum; 
        y_cm=ym_sum/m_sum;

        xc_box=(x0_box+x1_box)*0.5;
        yc_box=(y0_box+y1_box)*0.5;

        xguess=x_cm-xc_box;
        yguess=y_cm-yc_box;

        for i=1:wbox
            for j=1:wbox
                k=(i-1)*wbox+j;
                xymerge(k)=0;
                zmerge(k)=grab(i,j);
            end
        end 
      
        % Gaussian Fitting part
        
        beta0=[xguess,yguess,50,psf_w0/q,min(grab(:))];
        [betafit,resid,J,COVB,mse] = nlinfit(xymerge,zmerge,@gaussian_merge3,beta0);
        ci = nlparci(betafit,resid,'covar',COVB); 
        ci_err=(ci(:,2)-ci(:,1))/2; 
      

        yf=betafit(2)+yc_box;
        xf=betafit(1)+xc_box;
        a0=betafit(3);
        r0=abs(betafit(4));
        off=betafit(5);
        
        failed=0; 
        if(a0 < 0)
            n_fail_a0=n_fail_a0+1;
            failed=1;
        end
        if(xf > x1_box || xf < x0_box || yf > y1_box || yf < y0_box) 
           
            n_fail_outbox=n_fail_outbox+1;
            failed=1;
        end
        
        

        if(failed==0) 
       
            total_molecules=total_molecules+1;
            xcm_all(total_molecules)=x_cm;
            ycm_all(total_molecules)=y_cm;
            xf_all(total_molecules)=xf;
            yf_all(total_molecules)=yf;
            a0_all(total_molecules)=a0;
            r0_all(total_molecules)=r0;
            off_all(total_molecules)=off;
            
            % Saving info of each detected molecule
            
            frame_num_all(total_molecules)=fileloop;
            xf_err_all(total_molecules)=ci_err(1);
            yf_err_all(total_molecules)=ci_err(2);
            a0_err_all(total_molecules)=ci_err(3);
            r0_err_all(total_molecules)=ci_err(4);
            off_err_all(total_molecules)=ci_err(5);
            grab_sum_all(total_molecules)=grab_sum;
            grab_mol(:,:,total_molecules)=grab;         
            n_boxes=n_boxes+1; 
            boxes_xy(n_boxes,1)=high_pixel_x;
            boxes_xy(n_boxes,2)=high_pixel_y;
            boxes_xy(n_boxes,3)=1;
            box_matrix_ytop(total_molecules,1)=y0_box;
            box_matrix_xtop(total_molecules,1)=x0_box;

            
        end
    end
    
   
    imshow(i1,[]), axis image, colormap gray  
   title(['Frame: ' num2str(fileloop) ', ' num2str(total_molecules) ' molecules']); 
    hold on 
    draw_boxes(n_boxes,boxes_xy,rbox); %draw boxes is an outside function. 
    hold off %unfreeze the output
    drawnow; 
    
end


xcm_all=xcm_all(1:total_molecules);
ycm_all=ycm_all(1:total_molecules);
xf_all=xf_all(1:total_molecules);
yf_all=yf_all(1:total_molecules);
a0_all=a0_all(1:total_molecules);
r0_all=r0_all(1:total_molecules);
off_all=off_all(1:total_molecules);
frame_num_all=frame_num_all(1:total_molecules);
xf_err_all=xf_err_all(1:total_molecules);
yf_err_all=yf_err_all(1:total_molecules);
a0_err_all=a0_err_all(1:total_molecules);
r0_err_all=r0_err_all(1:total_molecules);
off_err_all=off_err_all(1:total_molecules);
grab_sum_all=grab_sum_all(1:total_molecules);
box_matrix_ytop=box_matrix_ytop(1:total_molecules);
box_matrix_xtop=box_matrix_xtop(1:total_molecules);
condition_sorting=zeros(total_molecules,1);

npix_all=pi*(r0_all.^2);    % area of molecule in square pixels
N=npix_all.*a0_all; % number of photons for each molecule

%Calculating the parameters 
            lp2=((r0_all*q).^2+(q^2)/12)*1./N+8*pi*((r0_all*q).^4)*(bkgn^2)/(q^2)*1./(N.*N);
            lp=sqrt(lp2); %loc prec in um 
       

     
           
             G_arr = [N npix_all xcm_all ycm_all xf_all yf_all a0_all r0_all off_all frame_num_all xf_err_all yf_err_all a0_err_all r0_err_all off_err_all grab_sum_all box_matrix_xtop box_matrix_ytop lp condition_sorting];
    G_arr1=G_arr;%copy of the inital without condition
        
        % Conditions chosen based on N HISTOGRAM plot
            G_arr(G_arr(:,1) <40,1) = 0;            
            G_arr(G_arr(:,1) >5000,1) = 0;            

        % |Defining the condition|
            TF1=G_arr(:,1)==0 ;
            G_arr(TF1,:)=[] ;   %|Deleting the npix_all and lp values based on zero N values|
            TF1c=TF1';
           %grab_mol_copy=grab_mol;
            grab_mol(:,:,TF1c)=[] ;   %|Deleting the npix_all and lp values based on zero N values|
          
            
         
           
        %|Condition on the area of molecule npix_all (in square pixels)|
            G_arr(G_arr(:,2) <9,2) = 0;             %# Set values in column 2 less than 9 (3x3 window) to 0
            G_arr(G_arr(:,2) >100,2) = 0;            %# Set values in column 2 greater than 25 (5x5 window) to 0

        % |Defining the condition|
            TF2=G_arr(:,2)==0 ;
            G_arr(TF2,:)=[] ; 
            TF2c=TF2';
             grab_mol(:,:,TF2c)=[] ; %|Deleting the N and lp values based on zero npix_all values|
        %-------------------------------------------------------------------------
                 
            G_arr(G_arr(:,19) >0.12,19) = 0;  
            TF3=G_arr(:,19)==0 ;
            G_arr(TF3,:)=[] ; 
            TF3c=TF3';
             grab_mol(:,:,TF3c)=[] ; 
        % Converting back to INDIVIDUAL SINGLE ARRAYS
            N=G_arr(:,1);
            npix_all=G_arr(:,2);
            xcm_all=G_arr(:,3);
            ycm_all=G_arr(:,4);
            xf_all=G_arr(:,5);
            yf_all=G_arr(:,6);
             
            a0_all=G_arr(:,7);
            r0_all=G_arr(:,8);
            off_all=G_arr(:,9);
            frame_num_all=G_arr(:,10);

            r0_err_all=G_arr(:,14);
            off_err_all=G_arr(:,15);
            grab_sum_all=G_arr(:,16);
           
            xf_err_all=G_arr(:,11);
            yf_err_all=G_arr(:,12);
            a0_err_all=G_arr(:,13);
            
            box_matrix_xtop=G_arr(:,17);
            box_matrix_ytop=G_arr(:,18);
            lp=G_arr(:,19);
            condition_sorting=G_arr(:,20);
        
% saving the file
    out_file = [an_dir,base_name,'_',num2str(n_start),'-',num2str(n_end),'-',num2str(min_thresh),'.mat'];
    save(out_file); 
   time=toc/60 
  
