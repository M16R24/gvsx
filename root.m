clc; clear; close;
tic;
x = linspace(0.001, 0.999, 100)';
T = (298:1:950)';
n = size(T);
dT = (T(end) - T(1) )/n(1);

ga = 131;
gb = 293;
R = 8.314;
gl = zeros(100, 1, n(1));
gbcc = zeros(100, 1, n(1));
ghcp = zeros(100, 1, n(1));

# loop to compute gibbs energy for each species

for j = 1:n(1)
    for i = 1:1000

        gl(:, 1, j) =  ((x*ga + (1-x)*gb) + (R*T(j)*(x.*log(x) + (1-x).*log(1-x))) + (x.*(1-x).*(-14935+10.371*T(j).*(x - (1-x))) + ...
            (-1789 + 1.143*T(j).*(x - (1-x)).^2) + (6533 - 6.591*T(j).*(x - (1-x)).^3)));
    
        gbcc(:,1, j) = (x*ga + (1-x)*gb) + (R*T(j)*(x.*log(x) + (1-x).*log(1-x))) + (x.*(1-x).*(-18335+8.49*T(j).*(x - (1-x))) + ...
            (3481.*(x - (1-x)).^2) + (2658 - 0.114*T(j).*(x - (1-x)).^3));
    
        ghcp(:,1, j) = (x*ga + (1-x)*gb) + (R*T(j)*(x.*log(x) + (1-x).*log(1-x))) + (x.*(1-x).*(6856.*(x - (1-x))) + ...
            (4000.*(x - (1-x)).^2) + (4000.*(x - (1-x)).^3));

        
    end
end


y = [gl(:, 1, 100), gbcc(:, 1, 100), ghcp(:, 1, 100)];
graph(x, y, 950);

# function to draw the graph

function graph(x_values, y_values, temperature)
    name = num2str(temperature);
    figure("Name", name);
    plot(x_values, y_values(:, 1))
    hold on
    plot(x_values, y_values(:, 2))
    plot(x_values, y_values(:, 3))
    hold off
    legend('g_{liq}', 'g_{bcc}', 'g_{hcp}');
    title("G vs X,T= " + name + "K");
    xlabel ('x_{B} \rightarrow')
    ylabel ('g J/mol', 'Interpreter', 'latex')
end

