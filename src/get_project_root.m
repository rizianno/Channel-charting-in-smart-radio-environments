function root = get_project_root()
% Returns the root folder of the project (Multi_Slot_Fusion_tsne)

    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);

    % Assuming this file is inside /src/
    root = fileparts(thisDir);

end