% test of matrix multiplication
% check.m

fida = fopen('A.txt','r'); % change filename if necessary
ma = fscanf(fida,'%d\n',1);
ka = fscanf(fida,'%d\n',1);
A = transpose(reshape(fscanf(fida,'%f\n',ka*ma),ka,ma));
fclose(fida);
%image(C);
%colormap(gray(256));
%set(gca,'DataAspectRatio',[1 1 1]);
fidb = fopen('B.txt','r'); % change filename if necessary
mb = fscanf(fidb,'%d\n',1);
kb = fscanf(fidb,'%d\n',1);
B = transpose(reshape(fscanf(fidb,'%f\n',kb*mb),kb,mb));
fclose(fidb);
C=A*B;