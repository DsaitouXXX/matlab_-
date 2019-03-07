function run_test(im_input,mean_num)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%

im_input=imresize(im_input,[256,256]);
%imshow(im_input)
im1=spe1(im_input);%����
figure,imshow(im1);
im2=im_input;%��ͼ

[x,y]=size(im1);

mean_num_des=sort(im1(:),'descend');%����ֵ 
label_all=(im1<=mean_num);
label_all=~fill_label(label_all);

for i=1:x
    for j=1:y
        if (label_all(i,j)==0) %meannum
            im2(i,j,:)=[255,255,255];
        end
    end
end
im2=uint8(im2);

figure,imshow(im2);
[im3,output_id]=mean_shift(im_input,50,70,0);% after fms
figure,imshow(im3);

for i=1:256
    for j=1:256
        if (im2(i,j,:)==[255,255,255])
            output_id(i,j)=0;
        end
    end
end
%%


%%%%%%%%%%%%%%%svm
im4=im2double(im3);
x=size(im3,1); y=size(im3,2);
cnt=0; limit1=2000; limit2=2000;

for i=1:x
    for j=1:y
        output_reshape((i-1)*x+j,:)=double(im_input(i,j,:));
        if (im1(i,j)>=mean_num_des(2000) && limit1>0)
            cnt=cnt+1; limit1=limit1-1;
            train_data(cnt,:)=output_reshape((i-1)*x+j,:);
            train_label(cnt)=1;
        end
        if (limit2>0 &&im1(i,j)<mean_num_des(2000))
                if ( output_id(i,j)==0 && ( (i>1 && output_id(i-1,j)~=output_id(i,j) ) ||  (j>1 && output_id(i,j-1)~=output_id(i,j) )  ))
                    cnt=cnt+1; limit2=limit2-1;
                    train_data(cnt,:)=output_reshape((i-1)*x+j,:);
                    train_label(cnt)=0;
                end
        end
    end
end

%%
model = svmtrain( train_data,train_label','kernel_function','linear');
predict_label = svmclassify(model,output_reshape);
%%��̬ѧ����
re_predict_label=reshape(predict_label,[256 256]);
re_predict_label=(re_predict_label');

for i=1:x
    for j=1:y
        if (re_predict_label(i,j)==0)
            label_all(i,j)=0;
        end
    end
end
label_all=~label_all;
figure,subplot(311),imshow(label_all),title('predict');
se1=strel('disk',2);%�����Ǵ���һ���뾶Ϊ2��ƽ̹��Բ�̽ṹԪ��
A2=imerode(label_all,se1);
subplot(312),imshow(A2),title('ʹ�ýṹԭʼdisk(2)��ʴ���ͼ��');
A3=imopen(A2,se1),subplot(313),imshow(A3),title('������'); 
im5=im_input; im6=im_input;
for i=1:x
    for j=1:y
       %if (~(predict_label((i-1)*x+j)==1 && output_id(i,j)==frequent_id) )
       %if ( output_id(i,j)~=frequent_id ) 
       if ( label_all(i,j)==1)
           im6(i,j,:)=[255,255,255];
       end
       if ( A3(i,j)==1 )
            im5(i,j,:)=[255,255,255];
        end
    end
end
%%
im5=uint8(im5);im6=uint8(im6);
figure;
subplot(121),imshow(im6),title('Ԥ��ͼ��');
subplot(122),imshow(im5),title('��̬ѧ������ͼ��');

end

