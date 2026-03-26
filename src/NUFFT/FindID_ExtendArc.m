function I = FindID_ExtendArc(M, pos_start, pos_end, x)
% FindID_ExtendArc

% Find indices I of x such that the the corresponding unit circle points
%  z = exp(2 * pi * x) are included in the arc
%  {exp(2 * pi * 1i * t / M): pos_start - 1 / 2 <= t < pos_end + 1 / 2}.

arguments (Input)
    M (1, 1) double;
    pos_start (1, 1) double;
    pos_end (1, 1) double;
    x (:, 1) double;
end

arguments(Output)
    I (:, 1) double;
end

half_length = 0.5 / M;
x_start = pos_start / M;
x_end = pos_end / M;

if pos_start == 0
    I = [find(x < x_end + half_length); find(x >= 1 - half_length)];
else
    I = intersect(...
        find(x < x_end + half_length), ...
        find(x >= x_start - half_length));
end

end