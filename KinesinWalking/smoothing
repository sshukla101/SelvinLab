% Written by Saurabh Shukla, Selvin lab, UIUC. 
% Last modified March 22, 2017.

# Go to the transformed folder of the analysis folder and slect a file. 
i=3

FileIn=dir('*.txt');
fid=fopen(FileIn(i).name);
Input = textscan(fid,'%f%f','CommentStyle','##');
yInput = Input{1};
position=yInput(1:2:end);
x=yInput(2:2:end);


%Wavelet denoising 
dn=wden(position,'modwtsqtwolog','s','mln',2,'sym4');
dn2=wden(position,'modwtsqtwolog','s','mln',3,'sym4');
dn3=wden(position,'modwtsqtwolog','s','mln',4,'sym4');

figure;
plot(dn)


set(gca,'ytick',linspace(0,160,21))
grid on
