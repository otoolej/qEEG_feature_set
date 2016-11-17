%---------------------------------------------------------------------
% generate EEG-like data (coloured Gaussian noise)
%---------------------------------------------------------------------
data_st=gen_test_EEGdata(5*60,64,1);


%---------------------------------------------------------------------
% define feature set (or can define in neural_parameters.m):
%---------------------------------------------------------------------
feature_set={'spectral_relative_power','rEEG_SD', 'connectivity_BSI'};
	
	
%---------------------------------------------------------------------
% estimate features:
%---------------------------------------------------------------------
feat_st=generate_all_features(data_st,[],feature_set);


