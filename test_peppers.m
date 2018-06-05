clear all; close all;
dbstop if error;

% obtain image
I = double(imread('peppers.png'));
figure;imshow(uint8(I));title('origin image');

% define foreground
[x,y] = meshgrid(142:150,133:433);
foreground = [x(:).';y(:).'];
[x,y] = meshgrid(210:218,71:465);
foreground = [foreground,[x(:).';y(:).']];

% defind background
[x,y] = meshgrid(30:38,133:433);
background = [x(:).';y(:).'];
[x,y] = meshgrid(364:372,133:433);
background = [background,[x(:).';y(:).']];

Iplot = I;
for i1 = 1:length(foreground)
    Iplot(foreground(1,i1),foreground(2,i1),:)=[255,0,0];
end
for i1 = 1:length(background)
    Iplot(background(1,i1),background(2,i1),:)=[0,0,255];
end
figure;imshow(uint8(Iplot));title('fore & back ground');

tic;
mattedImage = BP_matte(I,foreground, background);
toc;

