omega = [0.4781 ,  2.8867  ,  7.0404  ,  7.6839  , 14.1587];
freq = omega/(2*pi);
i = 5;
fid = fopen(['trans_D_end_',int2str(i),'.txt']);
D = fscanf(fid, '%f %f %f %f \n');

modal_analysis

figure
hold on
plot(abs(D/norm(D)),'x')
plot(abs(displ(:,i)/norm(displ(:,i))), 'o')
hold off

diff = abs(abs(displ(:,i)/norm(displ(:,i)))-abs(D/norm(D)))/max((displ(:,i)/norm(displ(:,i))));
fprintf('Max relative diff = %0.5f \n', 100*max(diff))