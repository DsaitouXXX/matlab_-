function [img1] = yz(im1,T)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
I=im2double(im1);
T=graythresh(I);
J=im2bw(I,T);
figure;
subplot(121),imshow(I);
subplot(122),imshow(J);
img1=J;
end

