function [img1] = yz(im1,T)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
I=im2double(im1);
T=graythresh(I);
J=im2bw(I,T);
figure;
subplot(121),imshow(I);
subplot(122),imshow(J);
img1=J;
end

