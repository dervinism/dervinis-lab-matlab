function test_writeUnitsTableH5()
%% test_writeUnitsTableH5.m — Tests for the low-level /units writer
%
% Builds a tiny synthetic NWB (MatNWB, for test-fixture setup only — the
% function under test never touches MatNWB), then writes /units via
% writeUnitsTableH5 and checks the result byte-for-byte against what was
% put in, including the two trickiest pieces: HDF5 object references
% (electrode_group column, spike_times_index target attribute) and the 2-D
% waveform_mean/waveform_sd dataspace orientation. Finishes with a pynwb
% read-back as the acceptance gate (pynwb is strict where MatNWB is not —
% see writeUnitsTableH5.m's header comment).
%
% Prerequisites on the MATLAB path:
%   - MatNWB (generateCore done) — test fixture only
%   - python + pynwb — the final read-back test is skipped without it

    fprintf('\n=== Running test_writeUnitsTableH5 ===\n\n');
    passed = 0; failed = 0;

    if exist('NwbFile', 'class') ~= 8
        fprintf('  SKIP: all tests (MatNWB not available)\n');
        return
    end

    work = tempname; mkdir(work);
    cleanup = onCleanup(@() rmdir(work, 's')); %#ok<NASGU>
    nwbFile = fullfile(work, 'test_units.nwb');
    local_build_fixture(nwbFile);

    %% Test 1: two units — every field round-trips correctly
    try
        units = local_two_unit_data();
        writeUnitsTableH5(nwbFile, units);

        colnames = cellstr(h5readatt(nwbFile, '/units', 'colnames'));
        assert(isequal(colnames(:), {'cluster_id';'type';'peak_channel_index'; ...
            'peak_channel_id';'x';'y';'z';'area';'lead_id';'electrode_group'; ...
            'spike_times';'waveform_mean';'waveform_sd'}), 'colnames order/content');
        assert(~any(strcmp(colnames, 'spike_times_index')), ...
            'spike_times_index must not be listed in colnames');

        id = h5read(nwbFile, '/units/id');
        assert(isequal(id, int64([0;1])), 'id');

        cid = h5read(nwbFile, '/units/cluster_id');
        assert(isequal(cid, [5;9]), 'cluster_id');

        typ = cellfun(@(x) strtrim(char(x)), h5read(nwbFile, '/units/type'), 'UniformOutput', false);
        assert(isequal(typ, {'unit';'mua'}), 'type');

        st = h5read(nwbFile, '/units/spike_times');
        assert(isequal(st, [0.1;0.2;0.3;0.5;0.7]), 'spike_times concatenation/order');

        sti = h5read(nwbFile, '/units/spike_times_index');
        assert(isequal(double(sti), [2;5]), 'spike_times_index cumulative counts');
        fid0 = H5F.open(nwbFile, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
        dsId0 = H5D.open(fid0, '/units/spike_times_index');
        typeId0 = H5D.get_type(dsId0);
        assert(H5T.get_sign(typeId0) == 0, ...
            'spike_times_index must be UNSIGNED (matches MatNWB; pynwb-validate rejects signed)');
        H5T.close(typeId0); H5D.close(dsId0); H5F.close(fid0);

        wm = h5read(nwbFile, '/units/waveform_mean');
        assert(isequal(size(wm), [3 2]), 'waveform_mean shape should be [nSamples x nUnits]');
        assert(isequal(wm, [1 10; 2 20; 3 30]), 'waveform_mean values/orientation');

        ws = h5read(nwbFile, '/units/waveform_sd');
        assert(isequal(size(ws), [2 2]), 'waveform_sd shape should be [nSamples x nUnits]');

        fid = H5F.open(nwbFile, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
        dsId = H5D.open(fid, '/units/electrode_group');
        refs = H5D.read(dsId);
        name1 = H5R.get_name(dsId, 'H5R_OBJECT', refs(:,1));
        name2 = H5R.get_name(dsId, 'H5R_OBJECT', refs(:,2));
        H5D.close(dsId);
        idxId = H5D.open(fid, '/units/spike_times_index');
        targetAttr = H5A.open(idxId, 'target');
        targetRef = H5A.read(targetAttr);
        targetName = H5R.get_name(idxId, 'H5R_OBJECT', targetRef);
        H5A.close(targetAttr); H5D.close(idxId);
        H5F.close(fid);
        assert(strcmp(name1, '/general/extracellular_ephys/LB_group_1'), 'unit1 electrode_group ref');
        assert(strcmp(name2, '/general/extracellular_ephys/LTAMicro_group_1'), 'unit2 electrode_group ref');
        assert(strcmp(targetName, '/units/spike_times'), 'spike_times_index target ref');

        fprintf('  PASS: Test 1 — two units round-trip correctly\n'); passed = passed + 1;
    catch e
        fprintf('  FAIL: Test 1 — %s\n', e.message); failed = failed + 1;
    end

    %% Test 2: single unit — waveform stays 2-D (no MatNWB singleton-flatten bug)
    try
        units = local_two_unit_data();
        units1 = local_row(units, 1);
        writeUnitsTableH5(nwbFile, units1);
        info = h5info(nwbFile, '/units/waveform_mean');
        assert(numel(info.Dataspace.Size) >= 2 && info.Dataspace.Size(2) == 1, ...
            sprintf('want 2-D [nSamples x 1], got %s', mat2str(info.Dataspace.Size)));
        fprintf('  PASS: Test 2 — single-unit waveform stays 2-D\n'); passed = passed + 1;
    catch e
        fprintf('  FAIL: Test 2 — %s\n', e.message); failed = failed + 1;
    end

    %% Test 3: idempotent re-run (no leftover / duplicate-group errors)
    try
        units = local_two_unit_data();
        writeUnitsTableH5(nwbFile, units);
        writeUnitsTableH5(nwbFile, units);
        id = h5read(nwbFile, '/units/id');
        assert(isequal(id, int64([0;1])), 'id after re-run');
        fprintf('  PASS: Test 3 — re-running is idempotent\n'); passed = passed + 1;
    catch e
        fprintf('  FAIL: Test 3 — %s\n', e.message); failed = failed + 1;
    end

    %% Test 4: pynwb can read the file back (the strict acceptance gate)
    try
        units = local_two_unit_data();
        writeUnitsTableH5(nwbFile, units);
        [ok, ~] = system('python -c "import pynwb"');
        if ok ~= 0
            fprintf('  SKIP: Test 4 — python/pynwb not available\n');
        else
            script = ['import pynwb\n' ...
                'with pynwb.NWBHDF5IO(r"' strrep(nwbFile,'\','\\') '", "r") as io:\n' ...
                '    nwbfile = io.read()\n' ...
                '    u = nwbfile.units[:]\n' ...
                '    assert len(u) == 2, len(u)\n' ...
                '    assert list(u["cluster_id"]) == [5.0, 9.0]\n' ...
                '    st = u["spike_times"]\n' ...
                '    assert list(st.iloc[0]) == [0.1, 0.2], st.iloc[0]\n' ...
                '    assert list(st.iloc[1]) == [0.3, 0.5, 0.7], st.iloc[1]\n' ...
                'print("PYNWB_OK")'];
            scriptFile = fullfile(work, 'read_units.py');
            fid = fopen(scriptFile, 'w'); fprintf(fid, '%s', sprintf(script)); fclose(fid);
            [st_, out] = system(sprintf('python "%s"', scriptFile));
            assert(st_ == 0 && contains(out, 'PYNWB_OK'), sprintf('pynwb read failed:\n%s', out));
            fprintf('  PASS: Test 4 — pynwb reads the file back correctly\n'); passed = passed + 1;
        end
    catch e
        fprintf('  FAIL: Test 4 — %s\n', e.message); failed = failed + 1;
    end

    %% Summary
    fprintf('\n=== Results: %d passed, %d failed ===\n\n', passed, failed);
    if failed > 0
        error('test_writeUnitsTableH5:failed', '%d test(s) failed', failed);
    end
end


%% ─────────────────────────── builders ──────────────────────────────────

function local_build_fixture(nwbOut)
    nwb = NwbFile('session_description','t','identifier','write-units-h5-test', ...
        'session_start_time', datetime(2025,1,1,'TimeZone','UTC'));

    nwb.general_devices.set('LB', types.core.Device('description','LB','manufacturer','Dixi'));
    nwb.general_devices.set('LTAMicro', types.core.Device('description','LTAMicro','manufacturer','Dixi'));
    nwb.general_extracellular_ephys.set('LB_group_1', types.core.ElectrodeGroup( ...
        'description','LB group 1','location','n/a','device', types.untyped.SoftLink('/general/devices/LB')));
    nwb.general_extracellular_ephys.set('LTAMicro_group_1', types.core.ElectrodeGroup( ...
        'description','LTAMicro group 1','location','n/a','device', types.untyped.SoftLink('/general/devices/LTAMicro')));

    % Column set/shape matches annotatedSz2nwb/test_check_nwb.m's local_build,
    % which is already proven to satisfy the stricter NWB>=2.9.0
    % ElectrodesTable schema pynwb enforces (a plain 4-6 column table is not
    % enough - see read_nwb_electrodes_table.m's docstring for the same
    % schema-version gotcha).
    nrow = 2;
    grp = [types.untyped.ObjectView('/general/extracellular_ephys/LB_group_1'); ...
           types.untyped.ObjectView('/general/extracellular_ephys/LTAMicro_group_1')];
    ChName = {'LB1';'LTAMicro1'};
    x = [1;10]; y = [2;11]; z = [3;12];
    gname = {'LB_group_1';'LTAMicro_group_1'};
    sc = @(v) string(compose("%d", v(:)));
    T = table(ChName, sc(1:nrow), sc(1:nrow), sc(1:nrow), sc(repmat(1,nrow,1)), sc(repmat(1,nrow,1)), ...
        sc(ones(nrow,1)), sc(ones(nrow,1)), sc(1:nrow), sc(zeros(nrow,1)), sc(zeros(nrow,1)), ...
        repmat("0",nrow,1), repmat("",nrow,1), repmat("",nrow,1), string(gname(:)), ...
        x, y, z, repmat("17",nrow,1), repmat("Left_Hippocampus",nrow,1), nan(nrow,1), ...
        repmat("Left_Hippocampus",nrow,1), repmat("none",nrow,1), grp, string(ChName(:)), ...
        string(gname(:)), repmat("L",nrow,1), repmat("n/a",nrow,1), ...
        'VariableNames', {'ChName','NLX_chNumber','leadSiteNumber','leadNumber', ...
          'NSitesThisLead','groupSize','NGroupsThisLead','iGroupThisLead','SiteNumberWithinGroup', ...
          'isMicro','isTetrode','isBad','associatedMacroLeadName','associatedMacroLeadNumber', ...
          'leadName','x','y','z','Destrieux_label','Destrieux_label_text','imp','location', ...
          'filtering','group','label','group_name','hemisphere','seizure_zone'});
    T.Properties.VariableDescriptions = repmat({'no description provided'}, 1, width(T));
    nwb.general_extracellular_ephys_electrodes = util.table2nwb(T, 'electrodes');

    etr = types.hdmf_common.DynamicTableRegion('table', ...
        types.untyped.ObjectView('/general/extracellular_ephys/electrodes'), ...
        'description','all','data', (0:nrow-1)');
    es = types.core.ElectricalSeries('data', int16(randi(100, nrow, 50)), ...
        'data_conversion', 1e-6, 'data_unit', 'volts', ...
        'electrodes', etr, 'starting_time', 0, 'starting_time_rate', 1000);
    nwb.acquisition.set('TimeSeries_1000_Hz', es);

    if isfile(nwbOut), delete(nwbOut); end
    nwbExport(nwb, nwbOut);
end


function units = local_two_unit_data()
    units = struct();
    units.clusterId          = [5; 9];
    units.type                = {'unit'; 'mua'};
    units.peakChannelIndex    = [1; 2];
    units.peakChannelId       = {'LB1'; 'LTAMicro1'};
    units.x                   = [1; 10];
    units.y                   = [2; 11];
    units.z                   = [3; 12];
    units.area                = {'Left_Hippocampus'; 'Left_Amygdala'};
    units.leadId               = {'LB'; 'LTAMicro'};
    units.electrodeGroupName  = {'LB_group_1'; 'LTAMicro_group_1'};
    units.spikeTimes          = {[0.1; 0.2], [0.3; 0.5; 0.7]};
    units.waveformMean        = [1 10; 2 20; 3 30];   % [nSamples=3, nUnits=2]
    units.waveformSd          = [0.1 1; 0.2 2];        % [nSamples=2, nUnits=2]
end


function out = local_row(units, i)
    out = struct();
    fn = fieldnames(units);
    for k = 1:numel(fn)
        name = fn{k};
        v = units.(name);
        if strcmp(name, 'spikeTimes')
            out.(name) = v(i);
        elseif strcmp(name, 'waveformMean') || strcmp(name, 'waveformSd')
            out.(name) = v(:, i);
        elseif iscell(v)
            out.(name) = v(i);
        else
            out.(name) = v(i);
        end
    end
end
