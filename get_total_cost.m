function [J] = get_total_cost(x, u, ts, N, l_cost_run, l_cost_final)
J = 0;
for k = 1:1:N
   l = full(l_cost_run(x(:,k),u(:,k)));
   J = J + l*ts;
end
l_f = full(l_cost_final(x(:,k)));
J = J + l_f;
end

