function shrinkplot(G, T)
%SHRINKPLOT Plots the objective funtion for the shrinking method.
%
%   By C. E. Mower, 12/11/14.
%

alpha_axis = linspace(0,1);

lambda_min = shrink_obj_func(alpha_axis, G, T);

figure;

h = plot(alpha_axis,lambda_min,'-k');
hold on;
plot([0 1], [0 0], '-k');

axis([0 1 lambda_min(1,1)-0.5 1]);

set(h, 'linewidth',2);

xL = get(gca,'XLim');

line(xL,[0 0], 'Color','k');

axis([0 1 lambda_min(1,1)-0.1 lambda_min(100)+0.1])

xlabel('0 \leq \alpha \leq 1', 'FontSize', 20);
ylabel('Minimum eigenvalue', 'FontSize', 20);
end
% -------------------------------------------------------------------------
function f = shrink_obj_func(alpha, A, T)
%SHRINK_OBJ_FUNC	Objective function for the Shrinking Method.

n = length(alpha);
f = zeros(1, n);

for i = 1:n
    S = alpha(i)*T + (1 - alpha(i))*A;
    f(i) = min(eig(S));
end

end
% -------------------------------------------------------------------------