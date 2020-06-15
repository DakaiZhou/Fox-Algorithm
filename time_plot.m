x = [2:7].^2;
y1 = [4.580951e-02, 3.760958e-01, 9.315650e-01, 2.028446e+00, 4.446542e+00, 8.138711e+00];
y2 = [1.323605e-02, 1.792560e-01, 2.504880e-01, 4.677510e-01, 1.108203e+00, 2.149676e+00]; 
%figure1 = figure;
semilogy(x, y1, x, y2, '--')
title('Weak scaling time estimates')
xlabel('Number of processes p')
ylabel('Time elapsed seconds')
legend('With I/O', 'Without I/O');