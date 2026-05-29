function w_cm_s = sinking_speed_named(d_cm, law_name)
% sinking_speed_named
% Return sinking speed in cm/s for named laws.

cfg = SimulationConfig();
law = lower(string(law_name));

switch law
    case "current"
        w_cm_s = local_current_law(d_cm, cfg);
    case "kriest_8"
        w_cm_s = (66 .* (d_cm .^ 0.62) .* 100) ./ cfg.day_to_sec;
    case "kriest_9"
        w_cm_s = (132 .* (d_cm .^ 0.62) .* 100) ./ cfg.day_to_sec;
    case "siegel_2025"
        d_mm = d_cm .* 10.0;
        w_cm_s = (20.2 .* (d_mm .^ 0.67) .* 100) ./ cfg.day_to_sec;
    otherwise
        error('Unknown law: %s', char(law_name));
end

w_cm_s(~isfinite(w_cm_s)) = 0;
w_cm_s(w_cm_s < 0) = 0;
end

function w_cm_s = local_current_law(d_cm, cfg)
r_v = 0.5 .* d_cm;
setcon = KernelLibrary.currentSetcon(cfg);
r_i = KernelLibrary.conservativeToFractalRadius(r_v, cfg);
w_cm_s = setcon .* (r_v .^ 3) ./ max(r_i, realmin);
end

