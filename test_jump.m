clear all; close all;
dbstop if error;

% obtain image
I = double(imread('jump2.png'));
figure;imshow(uint8(I));title('origin image');

% define foreground
[x,y] = meshgrid(84:120,90:94);
foreground = [x(:).';y(:).'];
[x,y] = meshgrid(52:56,87:93);
foreground = [foreground,[x(:).';y(:).']];

% defind background
[x,y] = meshgrid(24:185,18:29);
background = [x(:).';y(:).'];
[x,y] = meshgrid(25:182,141:151);
background = [background,[x(:).';y(:).']];

Iplot = I;
for i1 = 1:length(foreground)
    Iplot(foreground(1,i1),foreground(2,i1),:)=[255,0,0];
end
for i1 = 1:length(background)
    Iplot(background(1,i1),background(2,i1),:)=[0,0,255];
end
figure;imshow(uint8(Iplot));title('fore & back ground');

% tic;
% mattedImage = BP_matte(I,foreground, background);
% toc;