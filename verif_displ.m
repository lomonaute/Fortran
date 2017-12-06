
fid1 = fopen('displ_static.txt', 'r');
fid2 = fopen('displ_trans_0p5.txt', 'r');

static = fscanf(fid1, '%f');
trans = fscanf(fid2, '%f');

ind = ( static ~= 0);
mean(abs((trans(ind)-static(ind))./static(ind)))

%% displaccement of the mean line

nA = 2:3:62;
nAx = nA*2-1;
nAy = nA*2;

% plot(1:21, static(nAx),'ko', 1:21, trans(nAx), 'kx')
figure
plot(1:21, static(nAy),'ko:', 1:21, trans(nAy), 'kx')
ylabel('Displacement [m]')
yyaxis right
plot(1:21, abs((static(nAy)-trans(nAy))), '*-')
ylabel('Difference [m]')
xlim([1,21])

xlabel(' Node nr [-]')
legend ('Static', ' Transient') 