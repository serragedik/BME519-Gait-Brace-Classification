%% Batch gait-cycle segmentation for all Rokoko CSV trials
% This script processes all subject folders under one root folder.
% Outputs:
% - One *_segments.mat file per input CSV.
% - A summary table with cycle counts, QC numbers, and flags.
%
% Set plotOptions.generateTrialPlots = true to save QC plots.

clear; clc; close all;

%% CONFIG
fs = 100;  % Rokoko sampling frequency [Hz]

rootFolder = "C:\Users\user\Downloads\Computer Aided\Materials and Project\Yeni klasÃ¶r";
if strlength(rootFolder) == 0 || ~isfolder(rootFolder)
    selectedFolder = uigetdir(pwd, "Select root folder containing subject folders");
    if isequal(selectedFolder, 0)
        error("No root folder selected.");
    end
    rootFolder = string(selectedFolder);
end

outputFolder = fullfile(rootFolder, "_processed_segments");
segmentsFolder = fullfile(outputFolder, "segments_mat");
plotFolder = fullfile(outputFolder, "qc_plots");

timeCol = "Timestamp";         % Set [] if no timestamp column exists
flexCol = "RightKnee_flexion";
footYCol = "RightFoot_position_y";

flexFilter.order = 4;
flexFilter.cutoffHz = 6;
flexInvertSign = true;

footFilter.order = 4;
footFilter.cutoffHz = 4;

detect.eventMode = "postPeakMinimum";
detect.minHsDistanceSec = 0.65;
detect.minSwingPeakProminenceFrac = 0.12;
detect.minSwingPeakHeightFrac = 0.45;
detect.landingDescentFrac = 0.08;
detect.stanceThresholdFrac = 0.12;
detect.minBelowThresholdHoldSec = 0.08;
detect.landingSearchMaxSec = 0.90;
detect.refineToStablePlateau = true;
detect.stableSlopeFrac = 0.015;
detect.stableHoldSec = 0.04;
detect.maxStableSearchSec = 0.25;
detect.snapToPlateauStart = true;
detect.plateauSearchSec = 0.18;
detect.plateauBandFrac = 0.12;
detect.fillMissingCycles = true;
detect.maxCycleDurationMedianFactor = 1.45;
detect.insertSearchMarginSec = 0.20;
detect.minCycleDurationSec = 0.60;
detect.maxCycleDurationSec = 2.00;

nNormSamples = 101;

qc.removeOutlierCycles = true;
qc.minDurationMedianFactor = 0.65;
qc.maxDurationMedianFactor = 1.35;
qc.minCorrToMedianFlexion = 0.75;

summaryThresholds.minKeptCycles = 10;
summaryThresholds.maxRejectRate = 0.30;
summaryThresholds.maxCycleDurationCv = 0.20;
summaryThresholds.minMedianCorr = 0.80;

saveOptions.saveFilteredSignals = false;
saveOptions.saveAllCycles = true;

plotOptions.generateTrialPlots = false;
plotOptions.plotFlaggedOnly = true;
plotOptions.maxPlots = 30;
plotOptions.showRejectedCycles = true;

%% FIND FILES
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
if ~isfolder(segmentsFolder)
    mkdir(segmentsFolder);
end
if plotOptions.generateTrialPlots && ~isfolder(plotFolder)
    mkdir(plotFolder);
end

allFiles = dir(fullfile(rootFolder, "**", "*.csv"));
allFiles = allFiles(~contains(string({allFiles.folder}), string(outputFolder), "IgnoreCase", true));

if isempty(allFiles)
    error("No CSV files found under: %s", rootFolder);
end

fprintf("Found %d CSV files.\n", numel(allFiles));

%% PROCESS FILES
rows = cell(numel(allFiles), 1);
nPlotsSaved = 0;

for i = 1:numel(allFiles)
    sourceFile = string(fullfile(allFiles(i).folder, allFiles(i).name));
    subjectId = string(getParentFolderName(allFiles(i).folder));
    condition = parseConditionFromFileName(string(allFiles(i).name));

    fprintf("[%d/%d] %s | %s\n", i, numel(allFiles), subjectId, allFiles(i).name);

    try
        result = processOneTrial(sourceFile, subjectId, condition, fs, ...
            timeCol, flexCol, footYCol, flexFilter, footFilter, flexInvertSign, ...
            detect, qc, nNormSamples, saveOptions);

        resultFile = fullfile(segmentsFolder, subjectId, erase(string(allFiles(i).name), ".csv") + "_segments.mat");
        resultDir = fileparts(resultFile);
        if ~isfolder(resultDir)
            mkdir(resultDir);
        end

        save(resultFile, "result", "-v7.3");

        row = makeSummaryRow(result, sourceFile, resultFile, summaryThresholds, "");
        rows{i} = row;

        shouldPlot = plotOptions.generateTrialPlots && ...
            (~plotOptions.plotFlaggedOnly || row.flags ~= "OK") && ...
            nPlotsSaved < plotOptions.maxPlots;

        if shouldPlot
            plotFile = fullfile(plotFolder, subjectId + "_" + erase(string(allFiles(i).name), ".csv") + "_qc.png");
            saveTrialQcPlot(result, plotFile, plotOptions);
            nPlotsSaved = nPlotsSaved + 1;
        end

    catch ME
        rows{i} = makeErrorRow(sourceFile, subjectId, condition, ME.message);
        warning("Failed: %s\n%s", sourceFile, ME.message);
    end
end

summaryTable = struct2table(vertcat(rows{:}));
summaryCsv = fullfile(outputFolder, "segmentation_summary.csv");
summaryMat = fullfile(outputFolder, "segmentation_summary.mat");
writetable(summaryTable, summaryCsv);
save(summaryMat, "summaryTable");

fprintf("\nDone.\n");
fprintf("Segment files: %s\n", segmentsFolder);
fprintf("Summary CSV:   %s\n", summaryCsv);

disp("Flag counts:");
disp(groupsummary(summaryTable, "flags"));

%% LOCAL FUNCTIONS
function result = processOneTrial(sourceFile, subjectId, condition, fs, ...
    timeCol, flexCol, footYCol, flexFilter, footFilter, flexInvertSign, ...
    detect, qc, nNormSamples, saveOptions)

    T = readtable(sourceFile, "VariableNamingRule", "preserve");

    flexRaw = getNumericColumn(T, flexCol);
    footYRaw = getNumericColumn(T, footYCol);

    if flexInvertSign
        flexRaw = -flexRaw;
    end

    if isempty(timeCol)
        t = (0:numel(flexRaw)-1)' ./ fs;
    else
        tCandidate = getNumericColumn(T, timeCol);
        t = normalizeTimeVector(tCandidate, fs);
    end

    [bFlex, aFlex] = butter(flexFilter.order, flexFilter.cutoffHz/(fs/2), "low");
    flexFilt = filtfilt(bFlex, aFlex, flexRaw);

    [bFoot, aFoot] = butter(footFilter.order, footFilter.cutoffHz/(fs/2), "low");
    footYFilt = filtfilt(bFoot, aFoot, footYRaw);

    [eventIdxAll, swingPeakIdx, stanceThreshold] = detectLandingEventsFromFootY(footYFilt, fs, detect);

    eventIdx = removeImplausibleCycles(eventIdxAll, fs, ...
        detect.minCycleDurationSec, detect.maxCycleDurationSec);

    if detect.fillMissingCycles
        eventIdx = insertMissingCycleBoundariesFromFootY(footYFilt, eventIdx, fs, detect);
        eventIdx = removeImplausibleCycles(eventIdx, fs, ...
            detect.minCycleDurationSec, detect.maxCycleDurationSec);
    end

    if numel(eventIdx) < 2
        error("Fewer than 2 gait-cycle boundaries were detected.");
    end

    cycleDurSec = diff(eventIdx) ./ fs;
    [flexCyclesAll, footYCyclesAll, cycleTimePercent] = segmentAndNormalizeCycles( ...
        flexFilt, footYFilt, eventIdx, nNormSamples);

    if qc.removeOutlierCycles
        [goodCycleMask, cycleQc] = identifyGoodCycles(flexCyclesAll, cycleDurSec, qc);
    else
        goodCycleMask = true(1, size(flexCyclesAll, 2));
        cycleQc = emptyCycleQc(size(flexCyclesAll, 2));
    end

    result.sourceFile = sourceFile;
    result.subjectId = subjectId;
    result.condition = condition;
    result.fs = fs;
    result.nRows = height(T);
    result.trialDurationSec = (numel(flexRaw)-1) / fs;
    result.config.flexFilter = flexFilter;
    result.config.footFilter = footFilter;
    result.config.detect = detect;
    result.config.qc = qc;
    result.config.flexInvertSign = flexInvertSign;
    result.eventIdx = eventIdx;
    result.swingPeakIdx = swingPeakIdx;
    result.stanceThreshold = stanceThreshold;
    result.cycleDurSec = cycleDurSec;
    result.goodCycleMask = goodCycleMask;
    result.cycleQc = cycleQc;
    result.cycleTimePercent = cycleTimePercent;
    result.flexCycles = flexCyclesAll(:, goodCycleMask);
    result.footYCycles = footYCyclesAll(:, goodCycleMask);

    if saveOptions.saveAllCycles
        result.flexCyclesAll = flexCyclesAll;
        result.footYCyclesAll = footYCyclesAll;
    end

    if saveOptions.saveFilteredSignals
        result.t = t;
        result.flexRaw = flexRaw;
        result.flexFilt = flexFilt;
        result.footYRaw = footYRaw;
        result.footYFilt = footYFilt;
    else
        result.tStartEndSec = [t(1), t(end)];
    end
end

function row = makeSummaryRow(result, sourceFile, resultFile, thresholds, errorMessage)
    nCyclesRaw = numel(result.cycleDurSec);
    nCyclesKept = size(result.flexCycles, 2);
    nRejected = nCyclesRaw - nCyclesKept;

    if nCyclesRaw > 0
        rejectRate = nRejected / nCyclesRaw;
        meanDur = mean(result.cycleDurSec, "omitnan");
        stdDur = std(result.cycleDurSec, "omitnan");
        cvDur = stdDur / meanDur;
        minDur = min(result.cycleDurSec);
        maxDur = max(result.cycleDurSec);
    else
        rejectRate = nan;
        meanDur = nan;
        stdDur = nan;
        cvDur = nan;
        minDur = nan;
        maxDur = nan;
    end

    medianCorr = median(result.cycleQc.corrToTemplate(result.goodCycleMask), "omitnan");
    flags = makeFlags(nCyclesKept, rejectRate, cvDur, medianCorr, thresholds, errorMessage);

    row = baseSummaryRow(sourceFile, result.sourceFile, result.subjectId, result.condition);
    row.resultFile = string(resultFile);
    row.nRows = result.nRows;
    row.trialDurationSec = result.trialDurationSec;
    row.nEvents = numel(result.eventIdx);
    row.nCyclesRaw = nCyclesRaw;
    row.nCyclesKept = nCyclesKept;
    row.nRejected = nRejected;
    row.rejectRate = rejectRate;
    row.meanCycleDurSec = meanDur;
    row.stdCycleDurSec = stdDur;
    row.cvCycleDur = cvDur;
    row.minCycleDurSec = minDur;
    row.maxCycleDurSec = maxDur;
    row.medianCorrKept = medianCorr;
    row.flags = flags;
    row.errorMessage = string(errorMessage);
end

function row = makeErrorRow(sourceFile, subjectId, condition, errorMessage)
    row = baseSummaryRow(sourceFile, sourceFile, subjectId, condition);
    row.resultFile = "";
    row.nRows = nan;
    row.trialDurationSec = nan;
    row.nEvents = nan;
    row.nCyclesRaw = nan;
    row.nCyclesKept = nan;
    row.nRejected = nan;
    row.rejectRate = nan;
    row.meanCycleDurSec = nan;
    row.stdCycleDurSec = nan;
    row.cvCycleDur = nan;
    row.minCycleDurSec = nan;
    row.maxCycleDurSec = nan;
    row.medianCorrKept = nan;
    row.flags = "ERROR";
    row.errorMessage = string(errorMessage);
end

function row = baseSummaryRow(sourceFile, displayFile, subjectId, condition)
    row.sourceFile = string(sourceFile);
    row.displayFile = string(displayFile);
    row.subjectId = string(subjectId);
    row.fileName = string(getFileName(sourceFile));
    row.conditionCode = condition.conditionCode;
    row.condD = condition.condD;
    row.condE = condition.condE;
    row.condH = condition.condH;
    row.conditionSuffix = condition.suffix;
end

function flags = makeFlags(nCyclesKept, rejectRate, cvDur, medianCorr, thresholds, errorMessage)
    flagList = strings(0, 1);

    if strlength(string(errorMessage)) > 0
        flags = "ERROR";
        return;
    end

    if nCyclesKept < thresholds.minKeptCycles
        flagList(end+1) = "LOW_CYCLE_COUNT";
    end
    if rejectRate > thresholds.maxRejectRate
        flagList(end+1) = "HIGH_REJECT_RATE";
    end
    if cvDur > thresholds.maxCycleDurationCv
        flagList(end+1) = "HIGH_DURATION_VARIABILITY";
    end
    if medianCorr < thresholds.minMedianCorr
        flagList(end+1) = "LOW_MEDIAN_CORR";
    end

    if isempty(flagList)
        flags = "OK";
    else
        flags = strjoin(flagList, "|");
    end
end

function x = getNumericColumn(T, colName)
    names = string(T.Properties.VariableNames);
    colName = string(colName);
    idx = find(strcmpi(names, colName), 1);

    if isempty(idx)
        error("Column '%s' was not found.", colName);
    end

    x = T{:, idx};
    if iscell(x) || isstring(x) || ischar(x)
        x = str2double(string(x));
    end
    x = double(x(:));

    if all(isnan(x))
        error("Column '%s' could not be converted to numeric values.", colName);
    end

    x = fillmissing(x, "linear", "EndValues", "nearest");
end

function t = normalizeTimeVector(tCandidate, fs)
    tCandidate = double(tCandidate(:));
    tCandidate = fillmissing(tCandidate, "linear", "EndValues", "nearest");

    dtMed = median(diff(tCandidate), "omitnan");
    if dtMed > 0.5
        t = (tCandidate - tCandidate(1)) ./ 1000;
    elseif dtMed > 0
        t = tCandidate - tCandidate(1);
    else
        t = (0:numel(tCandidate)-1)' ./ fs;
    end
end

function [eventIdx, swingPeakIdx, stanceThreshold] = detectLandingEventsFromFootY(footY, fs, cfg)
    footY = footY(:);
    footRange = max(footY, [], "omitnan") - min(footY, [], "omitnan");
    low = prctile(footY, 10);
    high = prctile(footY, 90);

    stanceThreshold = low + cfg.stanceThresholdFrac * (high - low);
    minPeakProminence = cfg.minSwingPeakProminenceFrac * footRange;
    minPeakDistance = round(cfg.minHsDistanceSec * fs);
    minPeakHeight = low + cfg.minSwingPeakHeightFrac * (high - low);

    [~, swingPeakIdx] = findpeaks(footY, ...
        "MinPeakDistance", minPeakDistance, ...
        "MinPeakProminence", minPeakProminence, ...
        "MinPeakHeight", minPeakHeight);

    if isfield(cfg, "eventMode") && strcmpi(string(cfg.eventMode), "swingPeak")
        eventIdx = swingPeakIdx(:);
        return;
    end

    holdSamples = max(1, round(cfg.minBelowThresholdHoldSec * fs));
    landingSearchMaxSamples = round(cfg.landingSearchMaxSec * fs);
    eventCandidates = nan(numel(swingPeakIdx), 1);

    for i = 1:numel(swingPeakIdx)
        searchStart = min(swingPeakIdx(i) + 1, numel(footY));
        searchEnd = min(swingPeakIdx(i) + landingSearchMaxSamples, numel(footY));

        if i < numel(swingPeakIdx)
            searchEnd = min(searchEnd, swingPeakIdx(i+1) - round(0.10 * fs));
        end

        if searchEnd <= searchStart
            continue;
        end

        if isfield(cfg, "eventMode") && strcmpi(string(cfg.eventMode), "postPeakMinimum")
            candidate = detectPostPeakMinimum(footY, searchStart, searchEnd);
        elseif isfield(cfg, "eventMode") && strcmpi(string(cfg.eventMode), "fractionalDescent")
            candidate = detectLandingByFractionalDescent(footY, swingPeakIdx(i), searchStart, searchEnd, cfg);
        else
            below = footY(searchStart:searchEnd) <= stanceThreshold;
            firstSustainedBelow = findFirstRun(below, holdSamples);

            if isnan(firstSustainedBelow)
                candidate = nan;
            else
                candidate = searchStart + firstSustainedBelow - 1;

                if cfg.refineToStablePlateau
                    candidate = refineLandingToStablePlateau( ...
                        footY, candidate, searchEnd, fs, footRange, stanceThreshold, cfg);
                end

                if cfg.snapToPlateauStart
                    candidate = snapLandingToPlateauStart( ...
                        footY, candidate, searchEnd, fs, footRange, cfg);
                end
            end
        end

        if ~isnan(candidate)
            eventCandidates(i) = candidate;
        end
    end

    eventIdx = unique(eventCandidates(~isnan(eventCandidates)));
    eventIdx = eventIdx(:);
end

function landingIdx = detectPostPeakMinimum(footY, searchStart, searchEnd)
    landingIdx = nan;

    if searchEnd <= searchStart
        return;
    end

    idx = searchStart:searchEnd;
    [~, valleyRel] = min(footY(idx), [], "omitnan");
    landingIdx = idx(valleyRel);
end

function landingIdx = detectLandingByFractionalDescent(footY, peakIdx, searchStart, searchEnd, cfg)
    landingIdx = nan;

    if searchEnd <= searchStart
        return;
    end

    idx = searchStart:searchEnd;
    segment = footY(idx);
    [valleyValue, valleyRel] = min(segment, [], "omitnan");
    valleyIdx = idx(valleyRel);
    peakValue = footY(peakIdx);
    descentAmplitude = peakValue - valleyValue;

    if descentAmplitude <= 0
        return;
    end

    landingThreshold = valleyValue + cfg.landingDescentFrac * descentAmplitude;
    descentIdx = searchStart:valleyIdx;
    firstBelow = find(footY(descentIdx) <= landingThreshold, 1, "first");

    if ~isempty(firstBelow)
        landingIdx = descentIdx(firstBelow);
    end
end

function landingIdx = refineLandingToStablePlateau(footY, candidateIdx, searchEnd, fs, footRange, stanceThreshold, cfg)
    landingIdx = candidateIdx;
    refineEnd = min(searchEnd, candidateIdx + round(cfg.maxStableSearchSec * fs));

    if refineEnd <= candidateIdx
        return;
    end

    dy = [0; abs(diff(footY(:)))];
    stableSlopeThreshold = cfg.stableSlopeFrac * footRange;
    stableHoldSamples = max(1, round(cfg.stableHoldSec * fs));

    idx = candidateIdx:refineEnd;
    stableLow = footY(idx) <= stanceThreshold & dy(idx) <= stableSlopeThreshold;
    firstStableRun = findFirstRun(stableLow, stableHoldSamples);

    if ~isnan(firstStableRun)
        landingIdx = idx(firstStableRun);
    end
end

function landingIdx = snapLandingToPlateauStart(footY, candidateIdx, searchEnd, fs, footRange, cfg)
    landingIdx = candidateIdx;
    snapEnd = min(searchEnd, candidateIdx + round(cfg.plateauSearchSec * fs));

    if snapEnd <= candidateIdx
        return;
    end

    idx = candidateIdx:snapEnd;
    segment = footY(idx);
    localMin = min(segment, [], "omitnan");
    band = cfg.plateauBandFrac * footRange;
    firstNearPlateau = find(segment <= localMin + band, 1, "first");

    if ~isempty(firstNearPlateau)
        landingIdx = idx(firstNearPlateau);
    end
end

function firstIdx = findFirstRun(binaryVector, runLength)
    binaryVector = logical(binaryVector(:));
    firstIdx = nan;

    if numel(binaryVector) < runLength
        return;
    end

    runCount = conv(double(binaryVector), ones(runLength, 1), "valid");
    firstIdx = find(runCount == runLength, 1, "first");
end

function eventIdx = removeImplausibleCycles(eventIdxAll, fs, minDurSec, maxDurSec)
    eventIdxAll = eventIdxAll(:);
    if numel(eventIdxAll) < 2
        eventIdx = eventIdxAll;
        return;
    end

    keep = true(size(eventIdxAll));
    d = diff(eventIdxAll) ./ fs;

    shortGapIdx = find(d < minDurSec);
    keep(shortGapIdx + 1) = false;

    eventIdx = eventIdxAll(keep);

    if numel(eventIdx) >= 2
        d2 = diff(eventIdx) ./ fs;
        goodEvent = [true; d2 >= minDurSec & d2 <= maxDurSec];
        cycleGood = d2 >= minDurSec & d2 <= maxDurSec;
        eventUsed = false(size(eventIdx));
        eventUsed(1:end-1) = eventUsed(1:end-1) | cycleGood;
        eventUsed(2:end) = eventUsed(2:end) | cycleGood;
        eventIdx = eventIdx(eventUsed | goodEvent);
    end
end

function eventIdx = insertMissingCycleBoundariesFromFootY(footY, eventIdx, fs, cfg)
    eventIdx = sort(unique(eventIdx(:)));
    if numel(eventIdx) < 3
        return;
    end

    minDurSamples = round(cfg.minCycleDurationSec * fs);
    marginSamples = round(cfg.insertSearchMarginSec * fs);

    for iter = 1:4
        dSec = diff(eventIdx) ./ fs;
        plausibleDur = dSec >= cfg.minCycleDurationSec & dSec <= cfg.maxCycleDurationSec;
        if any(plausibleDur)
            medDur = median(dSec(plausibleDur), "omitnan");
        else
            medDur = median(dSec, "omitnan");
        end

        longGap = dSec > cfg.maxCycleDurationMedianFactor * medDur;
        if ~any(longGap)
            break;
        end

        newEvents = [];
        longGapIdx = find(longGap(:)');
        for k = longGapIdx
            left = eventIdx(k);
            right = eventIdx(k+1);
            searchStart = left + max(minDurSamples, marginSamples);
            searchEnd = right - max(minDurSamples, marginSamples);

            if searchEnd <= searchStart
                continue;
            end

            searchIdx = searchStart:searchEnd;
            if isfield(cfg, "eventMode") && strcmpi(string(cfg.eventMode), "swingPeak")
                [~, eventRel] = max(footY(searchIdx), [], "omitnan");
            else
                [~, eventRel] = min(footY(searchIdx), [], "omitnan");
            end
            candidate = searchIdx(eventRel);

            if candidate - left >= minDurSamples && right - candidate >= minDurSamples
                newEvents(end+1, 1) = candidate; %#ok<AGROW>
            end
        end

        if isempty(newEvents)
            break;
        end

        eventIdx = sort(unique([eventIdx; newEvents]));
    end
end

function [flexCycles, footYCycles, cyclePercent] = segmentAndNormalizeCycles(flex, footY, eventIdx, nNorm)
    nCycles = numel(eventIdx) - 1;
    flexCycles = nan(nNorm, nCycles);
    footYCycles = nan(nNorm, nCycles);
    cyclePercent = linspace(0, 100, nNorm)';

    for i = 1:nCycles
        idx1 = eventIdx(i);
        idx2 = eventIdx(i+1);
        idx = idx1:idx2;
        xOld = linspace(0, 100, numel(idx));

        flexCycles(:, i) = interp1(xOld, flex(idx), cyclePercent, "pchip");
        footYCycles(:, i) = interp1(xOld, footY(idx), cyclePercent, "pchip");
    end
end

function [goodMask, info] = identifyGoodCycles(flexCycles, cycleDurSec, qc)
    nCycles = size(flexCycles, 2);

    medDur = median(cycleDurSec, "omitnan");
    badDuration = cycleDurSec < qc.minDurationMedianFactor * medDur | ...
        cycleDurSec > qc.maxDurationMedianFactor * medDur;

    template = median(flexCycles, 2, "omitnan");
    corrToTemplate = nan(1, nCycles);

    for i = 1:nCycles
        x = flexCycles(:, i);
        valid = ~isnan(x) & ~isnan(template);
        if nnz(valid) < 10
            corrToTemplate(i) = nan;
        else
            c = corrcoef(x(valid), template(valid));
            corrToTemplate(i) = c(1, 2);
        end
    end

    badShape = corrToTemplate < qc.minCorrToMedianFlexion | isnan(corrToTemplate);
    goodMask = ~(badDuration(:)' | badShape);

    info.corrToTemplate = corrToTemplate;
    info.badDuration = badDuration(:)';
    info.badShape = badShape;
    info.nBadDuration = nnz(info.badDuration);
    info.nBadShape = nnz(info.badShape);
end

function info = emptyCycleQc(nCycles)
    info.corrToTemplate = nan(1, nCycles);
    info.badDuration = false(1, nCycles);
    info.badShape = false(1, nCycles);
    info.nBadDuration = 0;
    info.nBadShape = 0;
end

function condition = parseConditionFromFileName(fileName)
    [~, stem, ~] = fileparts(fileName);
    tokens = regexp(stem, "D(?<D>\d+)E(?<E>\d+)H(?<H>\d+)(?<suffix>.*)", "names");

    condition.conditionCode = string(stem);
    condition.condD = nan;
    condition.condE = nan;
    condition.condH = nan;
    condition.suffix = "";

    if ~isempty(tokens)
        condition.condD = str2double(tokens.D);
        condition.condE = str2double(tokens.E);
        condition.condH = str2double(tokens.H);
        condition.suffix = string(tokens.suffix);
    end
end

function folderName = getParentFolderName(folderPath)
    [~, folderName] = fileparts(folderPath);
end

function fileName = getFileName(filePath)
    [~, name, ext] = fileparts(filePath);
    fileName = name + ext;
end

function saveTrialQcPlot(result, plotFile, plotOptions)
    fig = figure("Visible", "off", "Color", "w", "Position", [100, 100, 1100, 850]);
    tiledlayout(1, 1, "TileSpacing", "compact", "Padding", "compact");
    nexttile;

    if plotOptions.showRejectedCycles && any(~result.goodCycleMask) && isfield(result, "flexCyclesAll")
        plot(result.cycleTimePercent, result.flexCyclesAll(:, ~result.goodCycleMask), ...
            "Color", [1.0 0.72 0.72], "LineWidth", 0.7); hold on;
    end

    plot(result.cycleTimePercent, result.flexCycles, "Color", [0.7 0.7 0.9], "LineWidth", 0.8); hold on;

    if ~isempty(result.flexCycles)
        plot(result.cycleTimePercent, mean(result.flexCycles, 2, "omitnan"), "b", "LineWidth", 2.2);
    end

    grid on;
    xlabel("Gait cycle (%)");
    ylabel("Right knee flexion");
    title(sprintf("%s | %s", result.subjectId, result.condition.conditionCode), "Interpreter", "none");
    exportgraphics(fig, plotFile, "Resolution", 160);
    close(fig);
end

