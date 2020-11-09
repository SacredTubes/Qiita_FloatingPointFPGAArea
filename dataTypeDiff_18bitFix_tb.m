%%
close all
clear
open_system('dataTypeDiff')
%%
bitWidth = 17;
scallingFactor = 3:(bitWidth-3);
% scallingFactor = bitWidth-3;  % 固定値でシミュレーション
dataPoint = 10000;

%% シミュレーションと誤差計算
for n = 1:numel(scallingFactor)
inDataMin(n) = -2^(bitWidth-1)/2^scallingFactor(n);
inDataMax(n) = (2^(bitWidth-1)-1)/2^scallingFactor(n);
powFactor(n) = nextpow2(inDataMax(n));
inData = linspace(inDataMin(n), inDataMax(n), dataPoint);
% Halfの最大値は[2^-14 2^16] (Doc記載）
% +-65504
% 正の最大値：half(2^16-17)。これ以上はInf扱い
% 負の最小値：half(-2^16+17)  同上

t = 0:length(inData)-1;
endTime = t(end);
simin0 = timeseries(inData, t);
% inDataMax = max(inData);
% inDataMin = min(inData);
%% シミュレーション実行
simDataMin = inDataMin(n);
simDataMax = inDataMax(n);
simBitWidth = bitWidth;

out = sim(gcs)

%% 差分を計算
% Simulinkデータを取得
outDouble = squeeze(out.logsout.getElement('out_double').Values.Data);
outSingle = squeeze(out.logsout.getElement('out_single').Values.Data);
outHalf = squeeze(out.logsout.getElement('out_half').Values.Data);
out18bit = squeeze(out.logsout.getElement('out_18bit').Values.Data);

% 差分
outDiffSingle = outDouble - double(outSingle);
outDiffHalf = outDouble - double(outHalf);
outDiff18bit = outDouble - double(out18bit);

% 標準偏差
stdSingle(n) = std(outDiffSingle)
stdHalf(n) = std(outDiffHalf)
std18bit(n) = std(outDiff18bit)

% 最大誤差
maxSingle(n) = max(abs(outDiffSingle))
maxHalf(n) = max(abs(outDiffHalf))
max18bit(n) = max(abs(outDiff18bit))

end

%% プロット
figure(1), hold on, grid on
set(gca, 'YScale', 'log')
plot(powFactor, stdSingle, 'm*',...
    powFactor, stdHalf, 'bx',...
    powFactor, std18bit, 'g^')
title('Standard Deviation')
xlabel('Data Range -2^n ~ 2^n')
ylabel('Std')
legend('Single', 'Half', '18bit','Location','northwest')


figure(2), hold on, grid on
set(gca, 'YScale', 'log')
plot(powFactor, maxSingle, 'm*',...
    powFactor, maxHalf, 'bx',...
    powFactor, max18bit, 'g^')
title('Max Error')
xlabel('Data Range -2^n ~ 2^n')
ylabel('Error')
legend('Single', 'Half', '18bit','Location','northwest')

%% 標準偏差をデータ範囲で正規化
stdSinglePersentage = stdSingle./(2.^powFactor)
stdHalfPersentage = stdHalf./(2.^powFactor)
std18bitPersentage = std18bit./(2.^powFactor)

meanStdSinglePersentage = mean(stdSinglePersentage)
meanStdHalfPersentage = mean(stdHalfPersentage)
meanStd18bitPersentage = mean(std18bitPersentage)

%% 最大偏差をデータ範囲で正規化
maxSinglePersentage = maxSingle./(2.^powFactor)
maxHalfPersentage = maxHalf./(2.^powFactor)
max18bitPersentage = max18bit./(2.^powFactor)

meanMaxSinglePersentage = mean(maxSinglePersentage)
meanMaxHalfPersentage = mean(maxHalfPersentage)
meanMax18bitPersentage = mean(max18bitPersentage)
