x =[ 0 0 1 0 1 1 0 1];

reading_data
%%
x = [ 0 1 1 0]';
y = [ 0 0 1 1]';
t_steps = length(P2);

%%
figure
hold on
hh1 = plot(x , y, 'kx:');
hh2 = plot(x , y, 'ko-');
xlim([-0.1, 1.4])
ylim([-0.1,1.2])
for t=1:20:t_steps/10
   set(hh1, 'XData',  x+5*P2(2*[1:4]-1,t)  , 'YData', y+5*P2(2*[1:4],t));
   %plot(x+10*P2(2*[1:4]-1,t)  , y+10*P2(2*[1:4],t), 'k.');
   drawnow
end
hold off