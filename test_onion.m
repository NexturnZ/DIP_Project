clear all; close all;
dbstop if error;

% obtain image
I = double(imread('onion.png'));
figure;imshow(uint8(I));title('origin image');

% define foreground
[x,y] = meshgrid(16:40,140:144);
foreground = [x(:).';y(:).'];


% defind background
[x,y] = meshgrid(96:100,15:184);
background = [x(:).';y(:).'];

Iplot = I;
for i1 = 1:length(foreground)
    Iplot(foreground(1,i1),foreground(2,i1),:)=[255,0,0];
end
for i1 = 1:length(background)
    Iplot(background(1,i1),background(2,i1),:)=[0,0,255];
end
figure;imshow(uint8(Iplot));title('fore & back ground');

tic;
mattedImage = BP_matte(I,foreground, background,2);
toc;
