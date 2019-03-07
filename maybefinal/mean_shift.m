function [output,output_id]=mean_shift(input_image,Hr,Th,sigma)

%计算中心点，算距离收敛
input1=input_image;
if (size(input1, 3) ~= 3)
    error('rgbhist:numberOfSamples', 'Input image must be RGB.')
end


input=imresize(input1,[256,256]);



color_bandwidth=Hr;
convergence_fact=Th;
%%%%%%%%%%%%%%%%% Color Histogram %%%%%%%%%%%%%%
fprintf('\n Computing Color histogram');
I=input; I_test=I;
nbins=256;
factor=256/nbins;

H1=zeros([nbins nbins nbins]);
ct=0;
for i=1:size(I,1)
    for j=1:size(I,2)
        
        p=double(reshape(I(i,j,:),[1 3])); %取某点位置的三色
        p=floor(p/factor)+1;%把颜色范围定义在nbins里
        ct=ct+1;%第ct个点
        coordinate(ct,:)=[i,j];
        arr1(ct,:)=[p(1),p(2),p(3)];

        H1(p(1),p(2),p(3))=H1(p(1),p(2),p(3))+1;  %标记颜色所属位置          
        pixels1(ct,:)=[I(i,j,:)];%第ct个点的三色位置    
    end
end

H1=H1(:);
H1=H1./sum(H1);
H_all=reshape(H1,[nbins,nbins,nbins]);

%%%%%%%%%%%%%%%%%%% Histogram done %%%%%%%%%%%%%%%%%%%%

ht=size(input,1);
wd=size(input,2);
output=input;
fprintf('\n Starting meanshift');
tic

mark_bins=zeros([nbins nbins nbins]);
center_sum=0;
%sigma=100;%分类的距离

for i=1:ct
    R=arr1(i,1); 
    G=arr1(i,2); 
    B=arr1(i,3);
    
    bin_r=floor(R/factor);
    bin_g=floor(G/factor); 
    bin_b=floor(B/factor);
    old_bins=[bin_r bin_g bin_b];  
    if (mark_bins(bin_r,bin_g,bin_b)>0) %1代表中心点 2代表被遍历过
        continue; 
    end

    %%%%%%%%% Meanshift part %%%%%%%%%%%
    
    dist=convergence_fact+1;
    color_bandwidth=floor(color_bandwidth/factor);
    
    while(dist>convergence_fact)
        
        hr=min(nbins,(bin_r+color_bandwidth)); 
        lr=max(1,(bin_r-color_bandwidth));
        
        hg=min(nbins,(bin_g+color_bandwidth)); 
        lg=max(1,(bin_g-color_bandwidth));
        
        hb=min(nbins,(bin_b+color_bandwidth));
        lb=max(1,(bin_b-color_bandwidth));
        
        s_r=0; s_b=0; s_g=0;
        weight=0;
    mark_bins(bin_r,bin_g,bin_b)=1;
    for k=lr:hr
        for l=lg:hg
            for m=lb:hb
                if (mark_bins(k,l,m)==0)
                    mark_bins(k,l,m)=2;
                     kk=norm((H_all(bin_r,bin_g,bin_b)-H_all(k,l,m))/color_bandwidth).^2;
                     gg=(1/sqrt(2*pi))*exp(-0.5*kk);
                    %gg=1;%高斯核函数 下次再挣扎吧（
                    s_r=s_r+k*H_all(k,l,m)*gg;
                    s_g=s_g+l*H_all(k,l,m)*gg;
                    s_b=s_b+m*H_all(k,l,m)*gg;
                    weight=weight+H_all(k,l,m)*gg;
                end
            end
        end
    end
        s_r=s_r/weight; s_g=s_g/weight; s_b=s_b/weight;
        rd=(s_r-bin_r); gd=(s_g-bin_g); bd=(s_b-bin_b);
        dist=sqrt(rd^2+gd^2+bd^2);
        bin_r=round(s_r); bin_g=round(s_g); bin_b=round(s_b);
    end
   %%%%%%%%%%%%%%%%%%% meanshift done %%%%%%%%%%%%%%%
   
   new_bins=[bin_r,bin_g,bin_b];
   xx=coordinate(i,1); yy=coordinate(i,2);
   input(xx,yy,:)=new_bins*factor;
	center_sum=center_sum+1;
	center_id(center_sum)=center_sum;
	center_map(center_sum)=i;
    center_mark(center_sum,:)=new_bins;
   if (center_sum>1)
        min_dist=sigma+1;
       for j=1:center_sum-1
           s_r=center_mark(j,1);
           s_g=center_mark(j,2);
           s_b=center_mark(j,3);
           rd=(s_r-bin_r); gd=(s_g-bin_g); bd=(s_b-bin_b);
           dist=sqrt(rd^2+gd^2+bd^2);
           if ( dist<sigma) 
			   if  (dist<min_dist)
			   		min_dist=sigma-dist;
			   		center_id(center_sum)=center_id(j);
               end
           end
       end
   end 
	
end


for i=1:center_sum
    fprintf('%d ',center_id(i));
end
fprintf('\n ');
for i=1:center_sum
    fprintf('%d ',center_map(i));
end
b_id=1:max(center_id);
c_id=histc(center_id,b_id);
[~,max_id_index]=max(c_id);
frequent_id=b_id(max_id_index);

fprintf('\ncenter_sum=%d\n frequent_id=%d\n',center_sum,frequent_id);
output_another=output;
for i=1:ht
    for j=1:wd
    
    p_y=i;
    p_x=j;
    
    R=input(p_y,p_x,1); G=input(p_y,p_x,2); B=input(p_y,p_x,3);
    
    bin_r=floor(double(R)/factor); 
    bin_g=floor(double(G)/factor); 
    bin_b=floor(double(B)/factor);
    old_bins=[bin_r bin_g bin_b];    
     
    min_dist=600; %sqrt 3*256^2
    k_now=1;
    for k=1:center_sum
    	s_r=center_mark(k,1); s_g=center_mark(k,2); s_b=center_mark(k,3);
        rd=(s_r-bin_r); gd=(s_g-bin_g); bd=(s_b-bin_b);
        dist=sqrt(rd^2+gd^2+bd^2);
        if (dist<min_dist || k==1)
        	min_dist=dist;
        	k_now=k;
        end
    end

    color_r=center_mark(center_id(k_now),1)*factor; %% computing color to be assigned
    color_g=center_mark(center_id(k_now),2)*factor;
    color_b=center_mark(center_id(k_now),3)*factor;
   
    output(i,j,1)=color_r;
    output(i,j,2)=color_g;
    output(i,j,3)=color_b;
    
    output_id(i,j)=center_id(k_now);
    if (center_id(k_now)==frequent_id)
      output_another(i,j,1)=0;
      output_another(i,j,2)=0;
      output_another(i,j,3)=0;
    else
        output_another(i,j,1)=255;
        output_another(i,j,2)=255;
        output_another(i,j,3)=255;
    end
    end    


end




