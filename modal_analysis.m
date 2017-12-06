modes = 5;
fid = fopen('modal_analysis.txt', 'r');

 A = fscanf(fid, '%f %f \n' );
 data = reshape(A, [1/modes*length(A),modes]);
 fclose all
 omega = data(2,:);
 displ = data(3:end, :);
 
%  % plot displ
%  nA = 2:3:62;
% nAx = nA*2-1;
% nAy = nA*2;
% figure
% plot(displ(nAx,:))
%  plot(displ(nAy,:), '*-')
%  xlabel(' Node nr [-]')
%  ylabel(' Vert. displacement [m]')
%  legend('Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Mode 5')
%  
% % Analytical
% fprintf('%0.6f',displ(:,1))