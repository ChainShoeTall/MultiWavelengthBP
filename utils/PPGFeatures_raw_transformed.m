classdef PPGFeatures
    methods (Static)

        function [flag1, flag2, peaks, valleys] = compute_cycle_pks_vlys(sig, fs, pk_th, removeStartEnd)
            if nargin < 4, removeStartEnd = true; end
            if nargin < 3, pk_th = 0.6; end

            % Detect peaks and valleys
            peaks, valleys = findpeaks_pyampd(sig, fs);
            % valleys = findpeaks_pyampd(max(sig)-sig, fs);

            % Clean start/end points
            peaks(peaks==1) = []; peaks(peaks==length(sig)) = [];
            valleys(valleys==1) = []; valleys(valleys==length(sig)) = [];

            if removeStartEnd
                if peaks(1) < valleys(1), peaks(1) = [];
                else, valleys(1) = []; end
                if peaks(end) > valleys(end), peaks(end) = [];
                else, valleys(end) = []; end
            end

            while ~isempty(peaks) && peaks(1) < valleys(1), peaks(1) = []; end
            while ~isempty(peaks) && peaks(end) > valleys(end), peaks(end) = []; end

            if isempty(peaks) || isempty(valleys)
                flag1 = true; flag2 = true; peaks = []; valleys = [];
                return;
            end

            meanVal = mean(sig(valleys));
            newPeaks = peaks(1);
            flag1 = false; flag2 = false;

            for i = 2:length(peaks)
                if (sig(peaks(i)) - meanVal) > (sig(newPeaks(end)) - meanVal)*pk_th
                    newPeaks(end+1) = peaks(i);
                end
            end
            if length(valleys)-1 ~= length(newPeaks), flag2 = true; end

            flag1 = ~isequal(peaks, newPeaks);
            peaks = newPeaks;
        end

        function [sp, dp, flag1, flag2, peaks, valleys] = compute_sp_dp(sig, fs, pk_th, removeStartEnd)
            [flag1, flag2, peaks, valleys] = ...
                PPGFeatures.compute_cycle_pks_vlys(sig, fs, pk_th, removeStartEnd);
            if ~isempty(peaks) && ~isempty(valleys)
                sp = median(sig(peaks));
                dp = median(sig(valleys));
            else
                sp = -1; dp = -1;
            end
        end

        function [cycles, peaksNorm, flag1, flag2, peaks, valleys] = extract_cycle_check(sig, fs, pk_th, removeStartEnd)
            [flag1, flag2, peaks, valleys] = ...
                PPGFeatures.compute_cycle_pks_vlys(sig, fs, pk_th, removeStartEnd);

            cycles = {};
            peaksNorm = [];
            if ~isempty(peaks) && ~isempty(valleys)
                for i = 1:(length(valleys)-1)
                    cycles{end+1} = sig(valleys(i):valleys(i+1));
                end
                if length(peaks)==length(valleys)-1
                    for i = 1:length(peaks)
                        peaksNorm(end+1) = peaks(i) - valleys(i);
                    end
                else
                    flag2 = true;
                end
            end
        end

        function [featNames, featAvg] = extract_feat_cycle(cycles, peaksNorm, fs)
            allFeats = [];
            for i = 1:length(cycles)
                try
                    [names, feats] = PPGFeatures.extract_temp_feat(cycles{i}, peaksNorm(i), fs);
                    featNames = names; allFeats(end+1,:) = feats; %#ok<AGROW>
                catch
                    warning('Ignoring faulty cycle.');
                end
            end
            if isempty(allFeats)
                featAvg = [];
            else
                featAvg = mean(allFeats,1);
            end
        end

        function [featNames, feat] = extract_feat_original(sig, fs, filtered, removeStartEnd)
            ppg = PPG_mat(sig, fs);
            [featuresMap, featNames, featCSV] = ppg.features_extractor(filtered, removeStartEnd);
            feat = cellfun(@(x)x(1), values(featuresMap));
        end

        function [freq, absFFT] = signal_fft(data, fs)
            L = length(data);
            Y = fft(data);
            P2 = abs(Y/L);
            freq = fs*(0:(L/2))/L;
            absFFT = P2(1:L/2+1);
        end

        function peaks = get_fft_peaks(absFFT, freq, minDist, numIter)
            if nargin<3, minDist=28; end
            if nargin<4, numIter=5; end
            locs = findpeaks_simple(absFFT(1:floor(end/6)), minDist);
            locs = locs(locs > minDist);
            peaks = locs(1:min(numIter,length(locs)));
        end

        function avgs = fft_peaks_neighbor_avg(absFFT, fftPeaks, win)
            if nargin<3, win=6; end
            avgs = arrayfun(@(p)mean(absFFT(max(1,p-win):min(end,p+win))), fftPeaks);
        end

        function cyclesStruct = extract_cycles_all_ppgs(waveforms, ppgPeaks, hrOffset, matchPosition, removeFirst)
            cyclesStruct = struct(); names = fieldnames(waveforms);
            offsets = round(hrOffset);
            if nargin<5, removeFirst=true; end

            if removeFirst, ppgPeaks = ppgPeaks(2:end-1); end
            for name=names'
                cyclesStruct.([name{1} '_cycles']) = {};
            end
            for p = ppgPeaks'
                lower = round(p - 0.25*offsets);
                upper = round(p + 0.75*offsets);
                if strcmp(matchPosition, 'dia_notches')
                    % skip for brevity
                end
                if lower<1 || upper>length(waveforms.(names{1})), continue; end
                for name=names'
                    cyclesStruct.([name{1} '_cycles']){end+1} = waveforms.(name{1})(lower:upper);
                end
            end
        end

        % -- Additional helper functions left out for brevity --
    end
end