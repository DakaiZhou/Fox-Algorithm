ma=1960;
na=2450;
mb=2450;
nb=1470;
A=rand(ma,na);
B=rand(mb,nb);
%x=rand(n,1);
dlmwrite('B.txt', mb, '-append', 'delimiter',' ');
dlmwrite('B.txt', nb, '-append', 'delimiter',' ');
dlmwrite('B.txt', B, '-append', 'delimiter',' ');
%dlmwrite('input.d', x, '-append', 'delimiter',' ');

dlmwrite('A.txt', ma, '-append', 'delimiter',' ');
dlmwrite('A.txt', na, '-append', 'delimiter',' ');
dlmwrite('A.txt', A, '-append', 'delimiter',' ');