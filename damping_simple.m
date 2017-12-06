% Textfile names

for i=1:9
    textfiles{i} = ['results_0p', int2str(i),'_0p', int2str(i),'.txt'];
end
textfiles{10} ='results_1p0_1p0.txt';

data = zeros(5,6000);

for i=1:10
 [T, data(i,:)] = reading_data(textfiles{i});
end
%%
dth = -0.358024369616;
figure
hold on
plot(T, data(2:2:10,:))
plot(T, dth*ones(size(T)),'k--')
hold off
ylabel('Displacement [m]')
xlabel('Time [s]')
legend(' \alpha =\beta =0.2', ' \alpha =\beta =0.4',' \alpha =\beta =0.6',' \alpha =\beta =0.8',' \alpha =\beta =1.0', 'Static')