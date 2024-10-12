% This program detects the fortunate molecules
% the output matlab file from molecule detection code should be given as the input for this
% code. So load that file in the workspace.

total_frames=max(frame_num_all);
total_molecules=size(frame_num_all,1);
molecules_per_frame=zeros(total_frames,1);
two_mol=zeros(total_molecules,2);
three_mol=zeros(total_molecules,3);
four_mol=zeros(total_molecules,4);
added_molecules=zeros(total_frames+3,1);

%FINDING ADDITONAL DATA WITH AVAILABLE DATA ---------------------------------------------------

for i=1:total_molecules
    mole=frame_num_all(i,1);
    molecules_per_frame(mole,1)= molecules_per_frame(mole,1)+1;
    
end

added_molecules(1,1)=molecules_per_frame(1,1);
for i=2:total_frames
    added_molecules(i,1)=molecules_per_frame(i,1)+added_molecules((i-1),1);
end

%% detection of repeated molecules based on centroid

f2=0;%f1,f2 are the variable which are used in IF condition to avoid recounting
f1=0;

for i=1:(total_molecules-molecules_per_frame(total_frames,1))
    z=frame_num_all(i,1);
    for j=added_molecules(z)+1:(added_molecules(z+1))
        difference_x =abs(xf_all(i)-xf_all(j));
        difference_y =abs(yf_all(i)-yf_all(j));
        if (sqrt(((difference_x).^2+(difference_y).^2)) <=r0_all(i,1))
            if(i~=f2)
            two_mol(i,1)=i;
            two_mol(i,2)=j;
            end
         end
     end     
end     
   



   two_a=two_mol(1:total_molecules,1);
   two_b=two_mol(1:total_molecules,2);
   
  
   two_a=nonzeros(two_a);
   two_b=nonzeros(two_b);
   
  
   num_two_repeated=size(two_a,1);
  
   final_two=zeros(num_two_repeated,2);
  
   
   final_two(1:num_two_repeated,1)=two_a;
   final_two(1:num_two_repeated,2)=two_b;
  
   f_2=two_b;
   f_1=two_a;
   
    for i=1:size(f_2,1)-1
        for j=i:size(f_2,1)-1
           di1=f_2(i,1)-f_2(j+1,1);
             di2=f_1(i,1)-f_1(j+1,1);
           if di1==0
               final_two(j,:)=0;
           end
           if di2==0
               final_two(j,:)=0;
           end
        end
    end 
   tf=final_two(:,1);
tf(tf>0)=1;
tf2=tf==0;
final_two(tf2,:)=[];

   untouch_final=final_two;
   fprintf(1,'Number of molecules repeated in two frames are %lu\n', num_two_repeated);
   
  
  %% The correlation part
 
cut_off=0.7;  % this is the correlation cutoff
final_two_copy=final_two;
final_corr=zeros(size(final_two_copy,1),3);
final_corr(:,1:2)=final_two_copy;
for i=1: size(final_two_copy,1)
    X1=1:1:512;
    Y1=X1;
    [Y,X]=meshgrid(X1,Y1);
    l=final_corr(i,1);
    k=final_corr(i,2);
    x01=G_arr(l,5);
    y01=G_arr(l,6);
    x02=G_arr(k,5);
    y02=G_arr(k,6);
    a01=G_arr(l,7);
    a02=G_arr(k,7);
    r01=G_arr(l,8);
    r02=G_arr(k,8);
    
    z1=a01*exp(-2*((X-x01).*(X-x01)+(Y-y01).*(Y-y01))/r01^2);
       z2=a02*exp(-2*((X-x02).*(X-x02)+(Y-y02).*(Y-y02))/r02^2);
       final_corr(i,3)=corr2(z1,z2);
end


final_corr_copy=final_corr;
tf=final_corr(:,3);

tf(tf<cut_off)=0;
tf(tf>cut_off)=1;

final_cor_mol=final_two_copy.*tf;

TF=final_cor_mol(:,1);
   TF(TF>0)=1;
   tf2=TF==0;
   final_cor_mol(tf2,:)=[];
   
%% Storing Molecule info for rendering the final image

  final_two=final_cor_mol;
   final_two_copy=final_two;
   g_sort=zeros(size(final_two,1),3);
   
for i=1 : size(final_two,1)
  a= final_two(i,2);
  if (a~=0)
  final_two(i,2)=0;
  [x1 y1]=find(final_two==a);
  final_two(i,2)=a;
    if (size(x1,1)~=0)
      k=3;
      c=x1(1,1);
      [final_two1,g_sort1]=repeated_molecules_sort(i,k,c,final_two,g_sort); % this is a function which sort molecules
      final_two=final_two1;
      g_sort=g_sort1;
    end
  end
end

%getting molecules from g sort
G_arr_copy=G_arr;
g_sort_copy=g_sort;
%final two
tf=g_sort(:,3);
tf(tf>0)=1;
tf2=tf==0;
final_two=final_two.*tf2;
tf=final_two(:,1);
tf(tf>0)=1;
tf2=tf==0;
final_two(tf2,:)=[];
%g_arr for 1 molecules


tf=zeros(size(G_arr,1),1);
for i=1: size(final_two_copy,1)
    k=final_two_copy(i,1);
    tf(k,1)=1;
    l=final_two_copy(i,2);
    tf(l,1)=1;
end
tf=not(tf);
  G_arr1=G_arr(:,:).*tf; 
  tf=G_arr1(:,1);
  tf(tf>0)=1;
  tf2=tf==0;
  G_arr1(tf2,:)=[];
  
  G_arr1(G_arr1(:,1) <90,1) = 0;             
    G_arr1(G_arr1(:,1) >5000,1) = 0;           
      TF1=G_arr1(:,1)==0 ;
      G_arr1(TF1,:)=[] ;  
      r0_all=G_arr1(:,8);
  N=G_arr1(:,1);
  lp2=((r0_all*q).^2+(q^2)/12)*1./N+8*pi*((r0_all*q).^4)*(bkgn^2)/(q^2)*1./(N.*N);
            lp=1.3*sqrt(lp2); %l
            G_arr1(:,19)=lp;
             G_arr1(G_arr1(:,19) >0.12,19) = 0;             
             
      TF1=G_arr1(:,19)==0 ;
      G_arr1(TF1,:)=[] ;  
            out_file = ['one_cen_',num2str(cut_off),'_',base_name,'.mat'];
    save(out_file,'G_arr1'); 
           
  
%G_arr for 2 molecule
G_arr1=G_arr;
for i=1:size(final_two)
    k=final_two(i,1);
    l=final_two(i,2);
    G_arr1(k,5)=(G_arr(k,3)+G_arr(l,3))/2; %xf average
    G_arr1(k,6)=(G_arr(k,4)+G_arr(l,4))/2; %yf average
    G_arr1(k,8)=(G_arr(k,8)+G_arr(l,8))/2; %r0 average
    G_arr1(k,1)=(G_arr(k,1)+G_arr(l,1)); %adding N
    G_arr1(k,20)=1; %condition 1 to preserve only this molecule for reconstruction
end

tf=G_arr1(:,20);
  tf(tf>0)=1;
  tf2=tf==0;
  G_arr1(tf2,:)=[];
  
  r0_all=G_arr1(:,8);
  N=G_arr1(:,1);
  lp2=((r0_all*q).^2+(q^2)/12)*1./N+8*pi*((r0_all*q).^4)*(bkgn^2)/(q^2)*1./(N.*N);
            lp=1.3*sqrt(lp2); %l
            G_arr1(:,19)=lp;
            out_file = ['two_cen_',num2str(cut_off),'_',base_name,'.mat'];
    save(out_file,'G_arr1'); 
%            
% g_arr for 3 molecules
g_sort(:,1:2)=final_two_copy;

tf=g_sort(:,3);
tf1=g_sort(:,4);
tf1(tf1>0)=1;
tf(tf>0)=1;
tf1=not(tf1);
tf=tf.*tf1;
tf2=tf==0;
g_sort_copy=g_sort;
g_sort_copy(tf2,:)=[];
final_three=g_sort_copy;
G_arr1=G_arr;
for i=1:size(final_three)
    k=final_three(i,1);
    l=final_three(i,2);
    m=final_three(i,3);
    G_arr1(k,3)=(G_arr(k,3)+G_arr(l,3)+G_arr(m,3))/3; %xcm average
    G_arr1(k,4)=(G_arr(k,4)+G_arr(l,4)+G_arr(m,4))/3; %ycm average
    G_arr1(k,8)=(G_arr(k,8)+G_arr(l,8)+G_arr(m,8))/3; %r0 average
    G_arr1(k,1)=G_arr(k,1)+G_arr(l,1)+G_arr(m,1); %adding N
    G_arr1(k,20)=1; %condition 1 to preserve only this molecule for reconstruction
end
tf=G_arr1(:,20);
  tf(tf>0)=1;
  tf2=tf==0;
  G_arr1(tf2,:)=[];
  
  r0_all=G_arr1(:,8);
  N=G_arr1(:,1);
  lp2=((r0_all*q).^2+(q^2)/12)*1./N+8*pi*((r0_all*q).^4)*(bkgn^2)/(q^2)*1./(N.*N);
            lp=1.3*sqrt(lp2); %l
            G_arr1(:,19)=lp;
            out_file = ['three_cen_',num2str(cut_off),'_',base_name,'.mat'];
    save(out_file,'G_arr1'); 
%            
            
 % G_arr for 4 molecules
    g_sort(:,1:2)=final_two_copy;

tf=g_sort(:,4);
tf(tf>0)=1;
tf1=g_sort(:,5);
tf1(tf1>0)=1;
tf1=not(tf1);
tf=tf1.*tf;
tf2=tf==0;
g_sort_copy=g_sort;
g_sort_copy(tf2,:)=[];
final_four=g_sort_copy;
G_arr1=G_arr;
for i=1:size(final_four)
    k=final_four(i,1);
    l=final_four(i,2);
    m=final_four(i,3);
    n=final_four(i,4);
    G_arr1(k,3)=(G_arr(k,3)+G_arr(l,3)+G_arr(m,3)+G_arr(n,3))/4; %xcm average
    G_arr1(k,4)=(G_arr(k,4)+G_arr(l,4)+G_arr(m,4)+G_arr(n,4))/4; %ycm average
    G_arr1(k,8)=(G_arr(k,8)+G_arr(l,8)+G_arr(m,8)+G_arr(n,8))/4; %r0 average
    G_arr1(k,1)=G_arr(k,1)+G_arr(l,1)+G_arr(m,1)+G_arr(n,1); %adding N
    G_arr1(k,20)=1; %condition 1 to preserve only this molecule for reconstruction
end
tf=G_arr1(:,20);
  tf(tf>0)=1;
  tf2=tf==0;
  G_arr1(tf2,:)=[];
  
  r0_all=G_arr1(:,8);
  N=G_arr1(:,1);
  lp2=((r0_all*q).^2+(q^2)/12)*1./N+8*pi*((r0_all*q).^4)*(bkgn^2)/(q^2)*1./(N.*N);
            lp=1.3*sqrt(lp2); %l
            G_arr1(:,19)=lp;  
            out_file = ['four_cen_',num2str(cut_off),'_',base_name,'.mat'];
    save(out_file,'G_arr1'); 
           
   % finding 5 and above molecules count
   tf=g_sort(:,5);
   tf(tf>0)=1;
   tf2=tf==0;
   g_sort(tf2,:)=[];
%    
    
