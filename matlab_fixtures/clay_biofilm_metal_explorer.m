% =========================================================================
% clay_biofilm_metal_explorer.m
% Read-Only MATLAB R2026b 4-Panel Clay-Biofilm-Metal Explorer (M7.9-CBM)
% Reads matlab_fixtures/cbm_m79_demo_bundle.json, validates schema/units,
% exports >= 1920x1080 PNG and JSON verification receipt.
% Strictly zero em-dashes and zero en-dashes.
% =========================================================================

function receipt = clay_biofilm_metal_explorer()
    fixture_dir = fileparts(mfilename('fullpath'));
    bundle_path = fullfile(fixture_dir, 'cbm_m79_demo_bundle.json');
    png_path    = fullfile(fixture_dir, 'clay_biofilm_metal_explorer_1920x1080.png');
    rcpt_path   = fullfile(fixture_dir, 'clay_biofilm_metal_explorer_receipt.json');

    val = validate_bundle_local(bundle_path);
    if ~val.valid
        error('SMBR:CBMBundleInvalid', 'Bundle validation failed: %s', val.error_code);
    end
    b = jsondecode(fileread(bundle_path));
    st = b.cbm_m79_stages;

    fig = figure('Visible', 'off', 'Color', [0.98 0.98 0.99], ...
                 'Position', [40, 40, 1920, 1080], 'InvertHardcopy', 'off');
    try
        theme(fig, 'light');
    catch
    end

    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle(tl, ...
        'M7.9 Clay-Biofilm-Metal Explorer (CBM-1) | Synthetic Demonstration | Read-Only 3D EDA', ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

    % Panel 1: Primary 5-Field State X = [theta_mobile, aw, Fe2, Fe3, S_aq] & d001 Hydration
    ax1 = nexttile(tl, 1);
    d001_nm = double(st.hydration_0W_to_3W_and_ectem.d001_nm(:));
    aw = double(b.fields.aw.inline_values(:));
    fe2 = double(b.fields.Fe2_mol_m3.inline_values(:));
    saq = double(b.fields.S_aq_mol_m3.inline_values(:));
    yyaxis(ax1, 'left');
    hB = bar(ax1, 1:numel(d001_nm), d001_nm, 0.55, 'FaceColor', [0.25 0.55 0.85], 'DisplayName', 'd001 (nm)');
    ylim(ax1, [0, 2.25]);
    ylabel(ax1, 'd001 basal spacing (nm)', 'Interpreter', 'none');
    yyaxis(ax1, 'right');
    hFe = plot(ax1, 1:numel(aw), fe2, '-o', 'LineWidth', 2.0, 'Color', [0.85 0.25 0.15], 'DisplayName', 'Fe2+ (mol/m^3)');
    hold(ax1, 'on');
    hSaq = plot(ax1, 1:numel(aw), saq * 5.0, '--s', 'LineWidth', 1.8, 'Color', [0.10 0.65 0.30], 'DisplayName', '5*S_aq (mol/m^3)');
    hold(ax1, 'off');
    ylim(ax1, [0, 0.78]);
    ylabel(ax1, 'Fe2+ (mol/m^3) and 5*S_aq', 'Interpreter', 'none');
    set(ax1, 'XTick', 1:5, 'XTickLabel', {'0W', '1W_Na', '2W_Ca', '3W_Na', 'EC-TEM'}, 'TickLabelInterpreter', 'none');
    grid(ax1, 'on');
    legend(ax1, [hB, hFe, hSaq], 'Location', 'northwest', 'Interpreter', 'none', 'FontSize', 8.5);
    title(ax1, sprintf('Panel 1 | Hydration 0W..3W & EC-TEM (Vac TEM: %s)', ...
          string(st.hydration_0W_to_3W_and_ectem.vacuum_tem_verdict)), 'Interpreter', 'none', 'FontSize', 10.5);

    % Panel 2: 5-Pool Sulfur Ledger & 4-Term Clay Density Separation
    ax2 = nexttile(tl, 2);
    pools_mM = double(st.five_pool_sulfur_ledger.pools_mol(:)) * 1e3;
    bar(ax2, 1:numel(pools_mM), pools_mM, 'FaceColor', [0.88 0.58 0.18]);
    set(ax2, 'XTick', 1:5, 'XTickLabel', {'SO4_aq', 'S_aq', 'FeS/Cu2S', 'Sorbed_S', 'H2S_gas'}, 'TickLabelInterpreter', 'none');
    grid(ax2, 'on');
    ylabel(ax2, 'Sulfur Inventory (mmol)');
    title(ax2, sprintf('Panel 2 | 5-Pool Sulfur Ledger (rel res = %.1e; Unresolved Density Triad: %s)', ...
          double(st.five_pool_sulfur_ledger.rel_residual), ...
          string(st.four_term_clay_density.unresolved_triad_verdict)), ...
          'Interpreter', 'none', 'FontSize', 10.5);

    % Panel 3: Discriminating Mechanism Gates CBM-G1, CBM-G2, CBM-G3
    ax3 = nexttile(tl, 3);
    axis(ax3, [0 1 0 1]);
    set(ax3, 'XTick', [], 'YTick', [], 'Box', 'on', 'Color', [0.98 0.98 0.99]);
    title(ax3, 'Panel 3 | Discriminating Gates CBM-G1, CBM-G2, CBM-G3', 'Interpreter', 'none', 'FontSize', 10.5);
    gTxt = sprintf([ ...
        'CBM-G1 (Abiotic Compaction / O2 Gate, Cerna 2026):\n' ...
        '  Sterile vs Biotic (1600 kg/m^3): %.1f vs %.1f um/yr\n' ...
        '  Abiotic dominated = %d | Biotic-MIC verdict: %s\n\n' ...
        'CBM-G2 (Screening vs Chelation vs Polymer Bridging):\n' ...
        '  Without empirical rheology : %s\n' ...
        '  With rheology + monovalent : %s\n\n' ...
        'CBM-G3 (Magnetite DIET vs Conductive e-Pili, Wang 2018):\n' ...
        '  Fe3O4 detected only        : %s\n' ...
        '  e-pili absent refutation   : %s'], ...
        double(st.discriminating_gates.CBM_G1_sterile_vs_biotic_um_yr(1)), ...
        double(st.discriminating_gates.CBM_G1_sterile_vs_biotic_um_yr(2)), ...
        logical(st.discriminating_gates.CBM_G1_abiotic_dominated), ...
        string(st.discriminating_gates.CBM_G1_biotic_mic_verdict), ...
        string(st.discriminating_gates.CBM_G2_without_rheology_verdict), ...
        string(st.discriminating_gates.CBM_G2_with_rheology_and_monovalent_verdict), ...
        string(st.discriminating_gates.CBM_G3_magnetite_only_verdict), ...
        string(st.discriminating_gates.CBM_G3_epili_absent_refutation_verdict));
    text(ax3, 0.5, 0.5, gTxt, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 9.0, 'FontName', 'FixedWidth', ...
         'BackgroundColor', [1.0 0.98 0.94], 'EdgeColor', [0.75 0.45 0.15], ...
         'Margin', 10, 'Interpreter', 'none');

    % Panel 4: Fisher Identifiability Spectrum & Bead-Calibrated FCM Observation Layer
    ax4 = nexttile(tl, 4);
    sv = double(st.identifiability_and_fcm.full_singular_values(:));
    semilogy(ax4, 1:numel(sv), max(sv, 1e-12), '-o', 'LineWidth', 2.2, 'Color', [0.15 0.50 0.80], ...
             'MarkerFaceColor', [0.15 0.50 0.80]);
    grid(ax4, 'on');
    xlim(ax4, [0.5, numel(sv) + 0.5]);
    xlabel(ax4, 'Fisher Information Singular Mode Index (Int64)');
    ylabel(ax4, 'Singular Value (log10 scale)');
    title(ax4, sprintf('Panel 4 | Fisher Rank %d/%d (Incomplete w/o aw+AVS: %d/%d) | FCM: %.2e cells/mL', ...
          int64(st.identifiability_and_fcm.full_fisher_rank), ...
          int64(st.identifiability_and_fcm.n_parameters), ...
          int64(st.identifiability_and_fcm.underdetermined_fisher_rank), ...
          int64(st.identifiability_and_fcm.n_parameters), ...
          double(st.identifiability_and_fcm.fcm_cells_per_mL)), ...
          'Interpreter', 'none', 'FontSize', 10.0);

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
    receipt.png_file = 'matlab_fixtures/clay_biofilm_metal_explorer_1920x1080.png';
    receipt.png_width_px = int64(infoImg.Width);
    receipt.png_height_px = int64(infoImg.Height);
    receipt.full_fisher_rank = int64(st.identifiability_and_fcm.full_fisher_rank);
    receipt.underdetermined_fisher_rank = int64(st.identifiability_and_fcm.underdetermined_fisher_rank);
    receipt.vacuum_tem_verdict = char(st.hydration_0W_to_3W_and_ectem.vacuum_tem_verdict);
    receipt.cbm_g1_verdict = char(st.discriminating_gates.CBM_G1_biotic_mic_verdict);
    receipt.cbm_g2_without_rheology_verdict = char(st.discriminating_gates.CBM_G2_without_rheology_verdict);
    receipt.cbm_g3_magnetite_only_verdict = char(st.discriminating_gates.CBM_G3_magnetite_only_verdict);
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
