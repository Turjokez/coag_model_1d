function out = size_class_general(Yfine, v_lower, diam_i, betas, N, edge_idx)
% size_class_general
% Group the fine sectional model into N coarse classes.
% Inputs are plain arrays only. No solver work is done here.
%
% Paper to code note:
% - Gelbard geometric sectionalization says many far-off gain terms are zero.
% - In this repo, b1 and b2 are gain terms into a fine section.
% - In this repo, b3, b4, and b5 are loss terms from a fine section.
% - This routine keeps the fine coag formulas, then just sums them into
%   coarse classes.

if nargin < 5
    error('size_class_general needs Yfine, v_lower, diam_i, betas, and N');
end
if nargin < 6
    edge_idx = [];
end

v_lower = v_lower(:);
diam_i = diam_i(:);
n_sec = numel(v_lower);

if n_sec < 2
    error('size_class_general needs at least 2 fine sections');
end
if numel(diam_i) ~= n_sec
    error('diam_i must have the same length as v_lower');
end
if any(~isfinite(v_lower)) || any(v_lower <= 0) || any(diff(v_lower) <= 0)
    error('v_lower must be positive and strictly increasing');
end
if ~isscalar(N) || ~isfinite(N) || N ~= round(N)
    error('N must be an integer scalar');
end
if N < 2 || N > n_sec
    error('N must be between 2 and the number of fine sections');
end
if ~isstruct(betas) && ~isa(betas, 'BetaMatrices')
    error('betas must be a struct or BetaMatrices object with b1..b5');
end
need_fields = {'b1','b2','b3','b4','b5'};
for i = 1:numel(need_fields)
    if ~isprop_or_field(betas, need_fields{i})
        error('betas is missing %s', need_fields{i});
    end
end

Yfine = normalize_state_matrix(Yfine, n_sec);
v_edges = [v_lower(1); 2.0 .* v_lower];
edge_idx = build_edge_idx(v_edges, N, edge_idx);
groups = [edge_idx(1:end-1)', edge_idx(2:end)' - 1];

[fine_gain, fine_loss] = compute_fine_coag_terms(Yfine, betas);
fine_net = fine_gain - fine_loss;

n_times = size(Yfine, 1);
C = zeros(n_times, N);
gain = zeros(n_times, N);
loss = zeros(n_times, N);

for ic = 1:N
    idx = groups(ic,1):groups(ic,2);
    C(:, ic) = sum(Yfine(:, idx), 2);
    gain(:, ic) = sum(fine_gain(:, idx), 2);
    loss(:, ic) = sum(fine_loss(:, idx), 2);
end

net = gain - loss;

diam_lo = diam_i(groups(:,1));
diam_hi = diam_i(groups(:,2));
diam_mid = sqrt(diam_lo .* diam_hi);

fine_total = sum(Yfine, 2);
coarse_total = sum(C, 2);
coarse_net_total = sum(net, 2);
fine_net_total = sum(fine_net, 2);

out = struct();
out.groups = groups;
out.edge_idx = edge_idx(:)';
out.C = C;
out.gain = gain;
out.loss = loss;
out.net = net;
out.diam_lo = diam_lo;
out.diam_hi = diam_hi;
out.diam_mid = diam_mid;
out.cons_mass_err = coarse_total - fine_total;
out.cons_rate_err = coarse_net_total - fine_net_total;
out.v_edges = v_edges;
out.v_edge_idx = v_edges(edge_idx);
out.fine_gain = fine_gain;
out.fine_loss = fine_loss;
out.fine_net = fine_net;
out.fine_total = fine_total;
out.coarse_total = coarse_total;
out.fine_net_total = fine_net_total;
out.coarse_net_total = coarse_net_total;
end

function tf = isprop_or_field(s, name)
if isstruct(s)
    tf = isfield(s, name);
else
    tf = isprop(s, name);
end
end

function val = get_prop_or_field(s, name)
if isstruct(s)
    val = s.(name);
else
    val = s.(name);
end
end

function Y = normalize_state_matrix(Yin, n_sec)
if isempty(Yin)
    error('Yfine cannot be empty');
end
if isvector(Yin)
    Yin = Yin(:)';
end
if size(Yin, 2) ~= n_sec
    error('Yfine must have %d fine sections in columns', n_sec);
end
if any(~isfinite(Yin(:)))
    error('Yfine must be finite');
end
Y = Yin;
end

function edge_idx = build_edge_idx(v_edges, N, edge_idx)
n_edges = numel(v_edges);
n_sec = n_edges - 1;

if ~isempty(edge_idx)
    edge_idx = edge_idx(:)';
    if numel(edge_idx) ~= N + 1
        error('explicit edge_idx must have length N+1');
    end
    if any(edge_idx ~= round(edge_idx))
        error('explicit edge_idx must be integer indices');
    end
    if edge_idx(1) ~= 1 || edge_idx(end) ~= n_edges
        error('explicit edge_idx must start at 1 and end at n_sections+1');
    end
    if any(diff(edge_idx) <= 0)
        error('explicit edge_idx must be strictly increasing');
    end
    if any(diff(edge_idx) < 1)
        error('each coarse class must contain at least one fine section');
    end
    return
end

log_edges = log10(v_edges);
target = linspace(log_edges(1), log_edges(end), N + 1);
edge_idx = ones(1, N + 1);
edge_idx(1) = 1;
edge_idx(end) = n_edges;

for k = 2:N
    [~, idx0] = min(abs(log_edges - target(k)));
    min_idx = k;
    max_idx = n_edges - (N + 1 - k);
    edge_idx(k) = min(max(idx0, min_idx), max_idx);
end

for k = 2:N
    edge_idx(k) = max(edge_idx(k), edge_idx(k-1) + 1);
end
for k = N:-1:2
    max_here = n_edges - (N + 1 - k);
    edge_idx(k) = min(edge_idx(k), max_here);
end
for k = 2:N
    if edge_idx(k) <= edge_idx(k-1)
        edge_idx(k) = edge_idx(k-1) + 1;
    end
end

if numel(unique(edge_idx)) ~= numel(edge_idx)
    error('automatic coarse edge build failed: repeated edges found');
end
end

function [gain, loss] = compute_fine_coag_terms(Yfine, betas)
n_times = size(Yfine, 1);
n_sec = size(Yfine, 2);

gain = zeros(n_times, n_sec);
loss = zeros(n_times, n_sec);

b1 = get_prop_or_field(betas, 'b1');
b2 = get_prop_or_field(betas, 'b2');
b3 = get_prop_or_field(betas, 'b3');
b4 = get_prop_or_field(betas, 'b4');
b5 = get_prop_or_field(betas, 'b5');

for it = 1:n_times
    v = Yfine(it, :)';
    v_pos = max(v, eps);
    v_r = v_pos';
    v_shift = [0, v_r(1:n_sec-1)];

    c_gain = v_r .* (v_r * b2) + (v_r * b1) .* v_shift;
    c_loss = v_r .* (v_r * (b3 + b4 + b5));

    gain(it, :) = c_gain;
    loss(it, :) = c_loss;
end
end
