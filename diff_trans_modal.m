omega = [0.4781 ,  2.8867  ,  7.0404  ,  7.6839  , 14.1587];
freq = omega/(2*pi);
i = 4;
[T, X] = reading_data(['trans_D_freq_',int2str(i),'.txt']);

%%
t = T(2)-T(1);             % Sampling period       
Fs = 1/t;            % Sampling frequency                    
L = length(T);             % Length of signal
%T = (0:L-1)*t;        % Time vector
n = 2^nextpow2(L);
Y = fft(X,n);

P2 = abs(Y/n);
P1 = P2(:,1:n/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);

f = 0: Fs/n: (Fs/2-Fs/n);
plot(f,P1(1:n/2)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 2])
P1(1:3)=0;
[a,b] = max(P1);
fprintf('Princ. frequency = %0.4f  +- %0.3f \n',f(b), Fs/n)
fprintf('Relative diff = %0.2f %% +- %0.2f %% \n',100*abs(f(b)-freq(i))/freq(i), 100*Fs/n/freq(i))
