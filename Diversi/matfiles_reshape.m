%% File to reshape mat files to avoid duplicating multiple times the variable pow_pkg

file_path = "../outputs/KalmanIdentRes/";
listing = dir(char(file_path));
file_names = "";
for indf = 1:length(listing)
    file_names = [file_names, string(listing(indf).name)];
end
files_data = file_path + file_names(contains(file_names,"KalmanIdentAll"))';
num_nodpkg = length(files_data);

for ind = 1:num_nodpkg
    data = load(files_data(ind));
    data.pow_pkg = data.pow_pkg{ind};
    save(files_data(ind), '-struct', 'data')
end
