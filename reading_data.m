%% This script reads data from a text file inside the Fortran folder

textfile = 'results.txt';

fid = fopen(textfile, 'r');
P=[];
T=[];
i=0;
tline = fgets(fid);
while tline~=-1 
if tline(1:3) == ' t='
    %disp('It works !')    
    C=textscan(tline, '%s %f');
    T= [T, C{2}];
else
    P = [P, cell2mat(textscan(tline, '%f %f %f %f'))];
    %P = [P, cell2mat(textscan(tline, '%f %f'))];
end
tline = fgets(fid);
end
fclose(fid);


%%
t_steps = length(T);
P2 = reshape(P, [length(P)/t_steps, t_steps]);

figure
plot(T, P2)
xlabel('Time [s]')
ylabel('Data analysed')

[~, t_bound] = min((sum(P2) - 10^50).^2);
fprintf('10^50 t_step = %0.0f \n', t_bound)