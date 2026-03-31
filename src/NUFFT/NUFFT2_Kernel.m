function k_z_w = NUFFT2_Kernel(z, w, N)

arguments (Input)
    z (:, 1) double;
    w (:, 1) double;
    N (1, 1) double;
end

arguments (Output)
    k_z_w (:, :) double;
end

phi = z ./ w.';
tmp_ind = abs(phi(:) - 1) <= eps;
k_z_w = (phi.^N - 1) ./ (phi - 1) / N;

% theta = angle(phi);
% k_z_w = expm1(theta * N * 1i) ./ expm1(theta * 1i) / N;

k_z_w(tmp_ind) = 1;

% DEBUG.
% m = size(z, 1);
% n = size(w, 1);
% k_z_w = zeros(m, n);
% for j = 1 : m
%     z_j  = z(j);
%     for k = 1 : n
%         w_k = w(k);
%         phi = z_j / w_k;
%         if abs(phi - 1) <= eps
%             v = 1;
%         else
%             v = (phi^N - 1) / (phi - 1) / N;
%         end
%         if abs(v) <= 1e-14
%             % figure();
%             % plot(real(z_j), imag(z_j), "Marker", "x");
%             % hold on
%             % plot(real(w_k), imag(w_k), "Marker", "+");
%             % keyboard;
%             v = 0;
%         end
%         k_z_w(j, k) = v;
%     end
% end

end
