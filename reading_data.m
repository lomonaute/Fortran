function [T, P2] = reading_data(textfile)

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
    %P = [P, cell2mat(textscan(tline, '%f %f %f %f'))];
    P = [P, cell2mat(textscan(tline, '%f'))];
end
tline = fgets(fid);
end
fclose(fid);


%%
t_steps = length(T);
P2 = reshape(P, [length(P)/t_steps, t_steps]);

end


