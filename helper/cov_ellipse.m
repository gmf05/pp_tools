function cov_ellipse(m_est,p)

N_cov = p.covariate_ind{end}(end);
% ind = 1:N_cov;
% 4-5 = noise
% 6-8 = laser
% 9-10 = interact
ind = [p.covariate_ind{3}(2) p.covariate_ind{4}(2)];
% ind = [p.covariate_ind{3} p.covariate_ind{4}];

% circle used to make ellipse
theta = 0:0.01:2*pi;
unit_circ = [cos(theta); sin(theta)];

% estimated parameters
b_hat = m_est.b(ind);
cov_mtx = m_est.W(ind,ind);
[e_vecs, e_vals] = eig(cov_mtx);
e_vals = sqrt(diag(e_vals));
P = size(cov_mtx,1);    %   # parameters

%   create scale+rotation matrix
A = zeros(P);
for p=1:P
    A(:,p) = 2*e_vals(p)*e_vecs(:,p);
end

%   ellipse = scaled + rotated circle
c_ellipse = A*unit_circ;
for p=1:P
    c_ellipse(p,:) = c_ellipse(p,:) + b_hat(p);
end

plot(b_hat(1), b_hat(2), 'rx', 'MarkerSize', 10, 'LineWidth',4); hold on
plot(c_ellipse(1,:), c_ellipse(2,:),'LineWidth',4);

% plot(exp(b_hat(1)), exp(b_hat(2)), 'rx', 'MarkerSize', 10, 'LineWidth',4); hold on
% plot(exp(c_ellipse(1,:)), exp(c_ellipse(2,:)),'LineWidth',4);

update_fig();

end