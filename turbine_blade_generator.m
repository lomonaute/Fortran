clc
clear

t_root = 2;  % blade tickness at the root
t_tip = 1;  % blade tickness at the tip
x0 = 0;  % distance from rotation axis to blade root
L_blade = 10;  % length of the blade

ne_x = 5;  % nodes in x-direction
ne_y = 5;  % nodes in y-direction
ne = ne_x * ne_y;

nn_x = ne_x  + 1;
nn_y = ne_y + 1;
nn = nn_x * nn_y;

x_vec = x0:L_blade/ne_x:L_blade+x0;
y_vec = t_root:-((t_root - t_tip)/ne_y):t_tip;
y_frac = 0:1/ne_y:1;

nodal_coordinates = zeros(nn, 3);
node_n = 1;
for i = 1 : ne_x + 1
    x = x_vec(i);
    for j = 1 : ne_y + 1
        y = y_frac(j) * y_vec(i);
        nodal_coordinates(node_n, :) = [node_n, x, y];
        node_n = node_n + 1;
    end
end
format_node = 'N, %i, %10.8f, %10.8f, 0. \n';
text_node = fprintf(format_node, nodal_coordinates');


fprintf('! Define element connectivity: EN, element #, nodal list \n')

% elements
element_coordinates = zeros(ne, 5);
element_n = 0;
for i = 1:ne_x
    for j = 1:ne_y
        n1 = (i - 1) * nn_y + j;
        n2 = i * nn_y + j;
        n3 = n2 + 1;
        n4 = n1 + 1;

        element_n = element_n + 1;
        element_coordinates(element_n, :) = [element_n, n1, n2, n3, n4]; 
    end
end

format_element = 'EN, %i, %i, %i, %i, %i \n';
text_element = fprintf(format_element, element_coordinates');

fprintf('! Define boundary support conditions: D, node #, dof label, value \n')

% boundary conditions
boundary_nodes = [(1:nn_y)', (1:nn_y)'];
format_boundary = 'D, %i, UX, 0.000000000\nD, %i, UY, 0.000000000 \n';
fprintf(format_boundary, boundary_nodes')

fprintf('! Define nodal load conditions: F, node #, dof label, value \n')
% input nodal loads

fprintf('! Define surface load conditions: SFE, element #, face #, PRES, 0, value \n')

el_num_lower = 1:ne_y:ne;
el_num_upper = el_num_lower + ne_y - 1;
el_num = [el_num_lower'; el_num_upper'];
%face
faces = [ones(ne_x, 1); 3 * ones(ne_x, 1)];
% value
pres_val = (x_vec(1:end-1) + x_vec(2:end))/2;  % value is equal to the length to the rotation axis
pres_val = [pres_val'; pres_val'];

format_pressure = 'SFE, %i, %i, %10.8f, 0, value \n';
fprintf(format_pressure, [el_num, faces, pres_val]')

fprintf('\n\nFINISH \n\n')




