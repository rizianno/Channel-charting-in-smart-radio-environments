function Data = load_dataset(params)
% LOAD_DATASET
% Loads the dataset specified in params.
%
% NOTE:
% Dataset field names may still use "RIS" (legacy naming).

    dataPath = params.dataPath;

    assert(isfile(dataPath), ...
        'Dataset not found at: %s\nCheck README.', dataPath);

    Data = load(dataPath);

end

