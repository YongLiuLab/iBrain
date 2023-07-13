function feature_mat = struct2mat(input_roi)

    for i = 1:length(input_roi.f1)
        if  input_roi.f1(i).maximum>0
            j1(i,1) = input_roi.f1(i).energy;
            j1(i,2) = input_roi.f1(i).kurtosis;
            j1(i,3) = input_roi.f1(i).maximum;
            j1(i,4) = input_roi.f1(i).mean;
            j1(i,5) = input_roi.f1(i).mad;
            j1(i,6) = input_roi.f1(i).median;
            j1(i,7) = input_roi.f1(i).minimum;
            j1(i,8) = input_roi.f1(i).range;
            j1(i,9) = input_roi.f1(i).rms;
            j1(i,10) = input_roi.f1(i).skewness;
            j1(i,11) = input_roi.f1(i).std;
            j1(i,12) = input_roi.f1(i).var;
            j1(i,13) = input_roi.f1(i).entropy;
            j1(i,14) = input_roi.f1(i).uniformity;

            j3(i,1)= input_roi.f3(i).Autocorrelation;
            j3(i,2)= input_roi.f3(i).ClusterProminence;
            j3(i,3)= input_roi.f3(i).ClusterShade;
            j3(i,4)= input_roi.f3(i).ClusterTendency;
            j3(i,5)= input_roi.f3(i).Contrast;
            j3(i,6)= input_roi.f3(i).Correlation;
            j3(i,7)= input_roi.f3(i).DifferenceEntropy;
            j3(i,8)= input_roi.f3(i).Dissimilarity;
            j3(i,9)= input_roi.f3(i).Energy;
            j3(i,10)= input_roi.f3(i).Entropy;
            j3(i,11)= input_roi.f3(i).Homogeneity1;
            j3(i,12)= input_roi.f3(i).Homogeneity2;
            j3(i,13)= input_roi.f3(i).IMC1;
            j3(i,14)= input_roi.f3(i).IMC2;
            j3(i,15)= input_roi.f3(i).IDMN;
            j3(i,16)= input_roi.f3(i).IDN;
            j3(i,17)= input_roi.f3(i).InverseVariance;
            j3(i,18)= input_roi.f3(i).MaximumProbability;
            j3(i,19)= input_roi.f3(i).SumAverage;
            j3(i,20)= input_roi.f3(i).SumEntropy;
            j3(i,21)= input_roi.f3(i).SumVariance;
            j3(i,22)= input_roi.f3(i).Variance;
            j3(i,23)= input_roi.f3(i).ShortRunEmphasis;
            j3(i,24)= input_roi.f3(i).LongRunEmphasis;
            j3(i,25)= input_roi.f3(i).GrayLevelNonuniformity;
            j3(i,26)= input_roi.f3(i).RunLengthNonuniformity;
            j3(i,27)= input_roi.f3(i).RunPercentage;
            j3(i,28)= input_roi.f3(i).LowGrayLevelRunEmphasis;
            j3(i,29)= input_roi.f3(i).HighGrayLevelRunEmphasis;
            j3(i,30)= input_roi.f3(i).ShortRunLowGrayLevelEmphasis;
            j3(i,31)= input_roi.f3(i).ShortRunHighGrayLevelEmphasis;
            j3(i,32)= input_roi.f3(i).LongRunLowGrayLevelEmphasis;
            j3(i,33)= input_roi.f3(i).LongRunHighGrayLevelEmphasis;
        
         end
    end
feature_mat = [j1,j3];
