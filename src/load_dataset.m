function Data = load_dataset(params)
% LOAD_DATASET
% Loads the dataset specified in params.
%
% NOTE:
% The dataset .mat file uses legacy field names (H_I_RIS1, H_o_RIS1, etc.).
% These are accessed directly here and mapped to EMS terminology in prepare_data.

    dataPath = params.dataPath;

    assert(isfile(dataPath), ...
        'Dataset not found at: %s\nCheck README.', dataPath);

    Data = load(dataPath);

end

