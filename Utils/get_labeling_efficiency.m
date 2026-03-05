function lab_eff = get_labeling_efficiency(use_simLabEff_modulation, phase_val,modmat_phase,modmat_lookup)
    if use_simLabEff_modulation
        phase_val = mod(phase_val, 2 * pi);
        phase_val = phase_val + (phase_val < -pi) * 2 * pi;
        phase_val = phase_val - (phase_val > pi) * 2 * pi;
        lab_eff = interp1(modmat_phase, modmat_lookup, phase_val, 'linear', 'extrap');
    else
        lab_eff = cos(phase_val);
    end
end