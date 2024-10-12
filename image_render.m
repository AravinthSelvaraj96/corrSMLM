% This program renders the final image. output from correlation part should be given as input

% GETTING THE SINGLE MOLECULE ARRAY DATA -------------------------------------
G_arr=G_arr1; %as the output of correlation part code saved as G_arr1
xw=size(iprod,2);%width of the image
yw=size(iprod,1);%height of the image


% INITIALIZING PARAMETERS ----------------------------------
expf=15;%expansion factor
q=0.06;%nm pixel size
A=1; % amplitude of gaussian function
cut_off=2; % cutoff in rendering molecule
image=zeros(yw*expf,xw*expf);
stop=size(G_arr,1);
box_size=7*expf; % 30x 30 gaussian will be generated
not_included=zeros(size(G_arr,1),1);
x=-box_size:1:box_size;
y=x;
[Y X]=meshgrid(x,y);
istart=1;
istop=stop;
% ---------------------------------------------------------

% GENERATING GAUSSIAN ERROR FUNCTION FOR SINGLE MOLECULES ----
for i=istart:istop
    xcm=round(G_arr(i,5)*expf);
    ycm=round(G_arr(i,6)*expf);
    lp=expf*(G_arr(i,19)/q);
    img_mol=A*exp(-((X.*X)+(Y.*Y))/(2*lp*lp));
    
    y1=ycm-box_size;
    y2=ycm+box_size;
    x1=xcm-box_size;
    x2=xcm+box_size;
    if(y1<1)
        y1=1;
    end
    if(y2>yw*expf)
        y2=yw*expf;
    end
    if(x1<1)
        x1=1;
    end
    if(x2>xw*expf)
        x2=xw*expf;
    end
    image(y1:y2,x1:x2)=image(y1:y2,x1:x2)+img_mol(1:y2-y1+1,1:x2-x1+1);
 end
% ---------------------------------------------------------


% SAVING & DISPLAYING THE RECONSTRUCTED IMAGE ------------- 
figure,imshow(image,[]);

pn=input('Enter file name : ','s');
pname=char(pn);
 outfile=strcat(pname,'_expf',num2str(expf),'_ar_A-',num2str(A),'_cut-',num2str(cut_off),'_mol,',num2str(istart),'-',num2str(istop),'.tif');
image1=(2^16)*(image/(max(max(image))));
  imwrite(uint16(image1),outfile,'Compression','none');
% ----------------------------------------------------------
