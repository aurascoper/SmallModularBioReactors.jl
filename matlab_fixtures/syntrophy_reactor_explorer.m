% =========================================================================
% syntrophy_reactor_explorer.m
% Read-Only MATLAB R2026b 4-Panel Syntrophy Reactor Explorer (SYN0..SYN6)
% Reads matlab_fixtures/syntrophy_demo_bundle.json, validates schema/units,
% exports >= 1920x1080 PNG and JSON verification receipt.
% Strictly zero em-dashes and zero en-dashes.
% =========================================================================

function receipt = syntrophy_reactor_explorer()
    fixture_dir = fileparts(mfilename('fullpath'));
    bundle_path = fullfile(fixture_dir, 'syntrophy_demo_bundle.json');
    png_path    = fullfile(fixture_dir, 'syntrophy_reactor_explorer_1920x1080.png');
    rcpt_path   = fullfile(fixture_dir, 'syntrophy_reactor_explorer_receipt.json');

    val = validate_bundle_local(bundle_path);
    if ~val.valid
        error('SMBR:SyntrophyBundleInvalid', 'Bundle validation failed: %s', val.error_code);
    end
    b = jsondecode(fileread(bundle_path));
    st = b.syntrophy_stages;

    fig = figure('Visible', 'off', 'Color', [0.98 0.98 0.99], ...
                 'Position', [40, 40, 1920, 1080], 'InvertHardcopy', 'off');
    try
        theme(fig, 'light');
    catch
    end

    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle(tl, ...
        'Syntrophy Reactor Explorer (SYN0 to SYN6) | Synthetic Demonstration Bundle | Read-Only View', ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

    % Panel 1: Partner Populations (qPCR X_donor vs X_methanogen & Methane Yield)
    ax1 = nexttile(tl, 1);
    Xd = double(b.fields.X_donor_kg_m3.inline_values(:));
    Xm = double(b.fields.X_methanogen_kg_m3.inline_values(:));
    CH4 = double(b.fields.CH4_produced_mol_m3.inline_values(:));
    x_idx = (1:numel(Xd))';
    bar(ax1, x_idx, [Xd, Xm, CH4], 'grouped');
    ylim(ax1, [0, 0.46]);
    grid(ax1, 'on');
    xlabel(ax1, 'Reactor Cell Index (Int64)');
    ylabel(ax1, 'Concentration (kg/m^3 or mol/m^3)');
    title(ax1, 'Panel 1 | SYN1 Partner-Specific qPCR Biomass & CH4 Yield (Total Biomass Refused)', ...
          'Interpreter', 'none', 'FontSize', 10.5);
    legend(ax1, {'X_donor (kg/m^3, qPCR)', 'X_methanogen (kg/m^3, qPCR)', 'CH4 produced (mol/m^3)'}, ...
           'Location', 'northeast', 'Interpreter', 'none', 'FontSize', 8.5);

    % Panel 2: Gas/Liquid Exchange & PHF-HUR Monod vs Moser Kinetics (70 C Cutoff)
    ax2 = nexttile(tl, 2);
    temps = double(st.SYN1_phf_hur_kinetics.temps_C(:));
    r_monod = double(st.SYN1_phf_hur_kinetics.monod_rates_mol_m3_s(:));
    r_moser = double(st.SYN1_phf_hur_kinetics.moser_rates_mol_m3_s(:));
    r_qmax = double(st.SYN1_phf_hur_kinetics.qmax_volumetric_mol_m3_s(:));
    hold(ax2, 'on');
    hQmax = plot(ax2, temps, r_qmax * 1e3, ':^', 'LineWidth', 1.8, 'Color', [0.40 0.40 0.40], ...
                 'DisplayName', 'PHF-HUR q_max capacity (peaks 65 C = 773 mmol/gVSS/d)');
    hM1 = plot(ax2, temps, r_moser * 1e3, '-o', 'LineWidth', 2.2, 'Color', [0.10 0.62 0.30], ...
               'DisplayName', 'Moser PHF-HUR at H2 = 0.02 mM (preferred fit)');
    hM2 = plot(ax2, temps, r_monod * 1e3, '--s', 'LineWidth', 2.0, 'Color', [0.20 0.45 0.85], ...
               'DisplayName', 'Monod baseline at H2 = 0.02 mM');
    xline(ax2, 70.0, 'r--', '70 C thermal cutoff (nil uptake)', 'LineWidth', 1.4, 'Interpreter', 'none');
    hold(ax2, 'off');
    grid(ax2, 'on');
    xlim(ax2, [25, 74]);
    xlabel(ax2, 'Temperature (deg C)');
    ylabel(ax2, 'Volumetric H2 Uptake (mmol / (m^3 * s))');
    title(ax2, sprintf('Panel 2 | SYN1-2 PHF-HUR Monod vs Moser & K_S Bias (gas/PHF-HUR = %.2fx)', ...
          double(st.SYN2_gas_liquid_ks_bias.ks_bias_ratio_gas_supplied)), 'Interpreter', 'none', 'FontSize', 10.5);
    legend(ax2, [hQmax, hM1, hM2], 'Location', 'northwest', 'Interpreter', 'none', 'FontSize', 8.0);

    % Panel 3: Spatial Syntrophy (5.30 -> 0.29 um Aggregate Fickian Diffusive Bound)
    ax3 = nexttile(tl, 3);
    d_um = double(st.SYN5_spatial_aggregate.distances_um(:));
    j_fick = double(st.SYN5_spatial_aggregate.fickian_flux_mol_m3_s(:));
    plot(ax3, d_um, j_fick, '-d', 'LineWidth', 2.2, 'Color', [0.45 0.25 0.78], ...
         'MarkerFaceColor', [0.45 0.25 0.78]);
    set(ax3, 'XDir', 'reverse');
    grid(ax3, 'on');
    xlabel(ax3, 'Calibrated Interspecies Distance (um, 5.30 to 0.29 um)');
    ylabel(ax3, 'Max Fickian H2 Diffusive Flux (mol / (m^3 * s))');
    title(ax3, sprintf('Panel 3 | SYN5 Coaggregation Diffusive Bound (%.2fx boost; CPM sites refused)', ...
          double(st.SYN5_spatial_aggregate.flux_amplification_5p30_to_0p29)), 'Interpreter', 'none', 'FontSize', 10.5);

    % Panel 4: Metabolic Exchange, Delta_r G', >390 ppm Butanol Shift & Evolutionary Quarantine
    ax4 = nexttile(tl, 4);
    axis(ax4, [0 1 0 1]);
    set(ax4, 'XTick', [], 'YTick', [], 'Box', 'on', 'Color', [0.98 0.98 0.99]);
    title(ax4, 'Panel 4 | SYN0-3-4 Ledger, Delta_r G'', Identifiability & Quarantine', ...
          'Interpreter', 'none', 'FontSize', 10.5);
    cardTxt = sprintf([ ...
        'SYN0 Stoichiometric Ledger:\n' ...
        '  Propionate e- rel residual : %.2e (closed = %d)\n' ...
        '  >390 ppm H2 butanol shift  : p_H2 = %.0f ppm (active = %d)\n\n' ...
        'SYN3 Thermodynamics & Regulatory Silence:\n' ...
        '  Delta_r G'' donor / meth   : %.1f / %.1f kJ/mol\n' ...
        '  Hamilton 2015 guard        : feasibility != monoculture flux\n\n' ...
        'SYN4 Identifiability & Evolutionary Quarantine:\n' ...
        '  Fisher rank (qPCR+PHF-HUR) : %d / %d (identifiable = %d)\n' ...
        '  Hydrogen Hypothesis        : %s (quarantined)\n' ...
        '  Asgard / Heimdallarchaeia  : %s (separate track)'], ...
        double(st.SYN0_stoichiometry.propionate_e_rel_residual), ...
        logical(st.SYN0_stoichiometry.propionate_e_closed), ...
        double(st.SYN0_stoichiometry.butanol_shift_pH2_ppm), ...
        logical(st.SYN0_stoichiometry.butanol_shift_active), ...
        double(st.SYN3_thermodynamics.delta_r_G_donor_J_mol) / 1000.0, ...
        double(st.SYN3_thermodynamics.delta_r_G_methanogen_J_mol) / 1000.0, ...
        int64(st.SYN4_identifiability_and_quarantine.fisher_rank), ...
        int64(st.SYN4_identifiability_and_quarantine.n_parameters), ...
        logical(st.SYN4_identifiability_and_quarantine.identifiable), ...
        string(st.SYN4_identifiability_and_quarantine.hydrogen_hypothesis_verdict), ...
        string(st.SYN4_identifiability_and_quarantine.asgard_heimdallarchaeia_verdict));
    text(ax4, 0.5, 0.5, cardTxt, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 9.0, 'FontName', 'FixedWidth', ...
         'BackgroundColor', [0.97 0.99 0.97], 'EdgeColor', [0.20 0.55 0.35], ...
         'Margin', 10, 'Interpreter', 'none');

    exportgraphics(tl, png_path, 'Resolution', 150);
    infoImg = imfinfo(png_path);
    if infoImg.Width < 1920 || infoImg.Height < 1080
        exportgraphics(tl, png_path, 'Resolution', 200);
        infoImg = imfinfo(png_path);
    end
    close(fig);

    receipt = struct();
    receipt.schema_version = '1.0.0';
    receipt.timestamp_utc = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
    receipt.bundle_valid = val.valid;
    receipt.png_file = 'matlab_fixtures/syntrophy_reactor_explorer_1920x1080.png';
    receipt.png_width_px = int64(infoImg.Width);
    receipt.png_height_px = int64(infoImg.Height);
    receipt.fisher_rank = int64(st.SYN4_identifiability_and_quarantine.fisher_rank);
    receipt.n_parameters = int64(st.SYN4_identifiability_and_quarantine.n_parameters);
    receipt.hydrogen_hypothesis_verdict = char(st.SYN4_identifiability_and_quarantine.hydrogen_hypothesis_verdict);
    receipt.asgard_heimdallarchaeia_verdict = char(st.SYN4_identifiability_and_quarantine.asgard_heimdallarchaeia_verdict);
    receipt.all_checks_passed = val.valid && (infoImg.Width >= 1920) && (infoImg.Height >= 1080);

    json_txt = jsonencode(receipt, 'PrettyPrint', true);
    fid = fopen(rcpt_path, 'w');
    assert(fid ~= -1, 'Could not write receipt');
    fwrite(fid, json_txt, 'char');
    fclose(fid);
end

function res = validate_bundle_local(bundle_fp)
    res = struct('valid', false, 'error_code', 'NONE');
    if ~isfile(bundle_fp)
        res.error_code = 'ERR_BUNDLE_FILE_NOT_FOUND';
        return;
    end
    b = jsondecode(fileread(bundle_fp));
    valid_units = {'1', 'occupancy_fraction', 'dimensionless', 'm', 'm/s', 'Pa', 'K', ...
                   'mol/m^3', 'mol/kgw', 'kg/m^3', 'Gy/s', 's^-1', 'J/mol', 'mol/(m^3*s)', 'mol/m^2'};
    fn = fieldnames(b.fields);
    for k = 1:numel(fn)
        fld = b.fields.(fn{k});
        if ~isfield(fld, 'unit') || ~any(strcmp(fld.unit, valid_units))
            res.error_code = 'ERR_INVALID_UNIT_CONTRACT';
            return;
        end
        if isfield(fld, 'inline_values') && any(~isfinite(double(fld.inline_values(:))))
            res.error_code = 'ERR_NAN_IN_NUMERIC_PAYLOAD';
            return;
        end
    end
    res.valid = true;
end
