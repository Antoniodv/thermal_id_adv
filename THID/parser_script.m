function res = parser_script(obj, ismetric, test_idx, tconf_idx, tround, dfmat_idx, dtype_idx, path_to_pm_file)
%%%%
if ((nargin < 5) || isempty(dfmat_idx) )
	dfmat_idx = 1;
end
if ((nargin < 6) || isempty(dtype_idx))
	dtype_idx = 1;
end

%% Select Tests
% BASE_DIR = "C:\Users\anton\OneDrive\Desktop\phd\LEONARDO\SEED\thermal_identification";
% DELIV_FOLDER_NAME = "DEC2_deliv";
% TEST_NAME = "DEC2_prbs_random";
% N_CORES = "4";
% N_ROUND = "0";

% SLASH = "\";
% 
% FOLDER_NAME = "prbs_112_no_norm"
% FILE_NAME = "FEB5_prbs_random.pkl"
% 
% 
% file_path = strcat(BASE_DIR, SLASH, "Results", SLASH, FOLDER_NAME, SLASH, FILE_NAME)

file_path = path_to_pm_file;


% disp(file_path)

% pyenv('Version', 'C:\Users\anton\AppData\Local\Programs\Python\Python311\python.exe')
% py.print('PYTHON CORRECTLY CALLED BY MATLAB');


%% Extract:

switch dfmat_idx
	case 1 % .pkl
		try
			fid=py.open(file_path,'rb');
		catch ME
			res = [];
            disp("parser_script: no file found at the pointed path")
			return;
		end
		data=py.pickle.load(fid);

		%{
		% Replace invalid characters in keys with underscores
		new_data = containers.Map;
		dict_keys = data.keys();
		for i = 1:length(dict_keys)
    		key = char(dict_keys{i});
    		valid_key = regexprep(key, '[^a-zA-Z0-9_]', '_');  % Replace invalid characters with _
    		new_data(valid_key) = data{dict_keys{i}};
		end
		%}

		switch dtype_idx
			case 1 %list of Dictionary
				% Convert keys and values to MATLAB-compatible cell arrays
				%keys = cell(py.list(data.keys()));
				%dict_keys = cellfun(@(x) string(char(x)), keys, 'UniformOutput', false); %cell(py.list(data.keys()));     % Convert keys to MATLAB cell array
				dict_keys = py.list(data.keys());
				dict_values = py.list(data.values());       % List of values from Python dictionary
				
				% Initialize cell array for converted values
				mat_data = cell(size(dict_values));
				
				% Convert each DataFrame value to a MATLAB matrix
				for i = 1:length(dict_values)
					python_df = dict_values{i};
					matlab_matrix = double(python_df.values);  % Convert Python DataFrame to MATLAB matrix
					mat_data{i} = matlab_matrix;
                end

				mat_data = mat_data';
				dict_keys = string(dict_keys)';

				% Combine keys and MATLAB matrices into a cell array
				%mat_data = [dict_keys', mat_data'];  % Transpose to create a 2-column cell array
                
			case 2 %Dictionary/Dataframe
				% Convert keys and values to MATLAB-compatible cell arrays
				%keys = cell(py.list(data.keys()));
				%dict_keys = cellfun(@(x) string(char(x)), keys, 'UniformOutput', false); %cell(py.list(data.keys()));     % Convert keys to MATLAB cell array
				% Assuming `data` is a Python Pandas DataFrame

				numpy_array = data.values;          % Extract the NumPy array from the DataFrame
                mat_data = double(numpy_array);     %double(py.array.array('d', numpy_array.flatten()));  % Convert to MATLAB double array
				
                %mat_data = reshape(mat_data, size(numpy_array));  % Reshape to match the original DataFrame dimensions
			case 3 %list
				mat_data = cell(data);
			case 4 %Numpy Array
				mat_data = double(data);

			otherwise
				error("Wrong Python Data Type selected!");
		end
	case 2 %.csv
		%readmatrix %TODO
	otherwise
		error("Wrong Data Format selected!");
end


%% ...

if ismetric==1
	res = [];
	
	cores_tot_num = size(mat_data{1},2);
	if mod(cores_tot_num, obj.num_socket) ~= 0
		error("Cores division into socket is not uniform, please change the code");
	end
	cores_num = cores_tot_num / obj.num_socket;
	
	temp = mat_data{dict_keys == obj.TEMP_CORE};
	
	%C0 = mat_data{dict_keys == APERF_CORE};
	mperf = mat_data{dict_keys == obj.MPERF_CORE};
	freqd = mat_data{dict_keys == obj.FREQ};
	clk = mat_data{dict_keys == obj.CLK_CORE};
	clk_ref = mat_data{dict_keys == obj.CLK_REF_CORE};
	tsc = mat_data{dict_keys == obj.TSC};
	top_down = mat_data{dict_keys == obj.TOP_DOWN};
	

	sz = size(mat_data{dict_keys == obj.C_CORE(1)});
	c = nan([sz, length(obj.C_CORE)]);

    
	%TODO: ........
	minsz = min(size(mat_data{dict_keys == obj.C_CORE(1)},1), size(tsc,1));
	tsc = tsc(1:minsz,:);
	for i=1:length(obj.C_CORE)
		capa = mat_data{dict_keys == obj.C_CORE(i)};
		capa = capa(1:minsz,:);
		c(:,:,i) = capa ./ tsc;		
	end
	for i=1:size(c,3)
		total_cstate = sum(c,i);
	end
	
	%%% CPU TODO CHANGE
	C2_unc = mat_data{dict_keys == "C2_cpu"};
	C3_unc = mat_data{dict_keys == "C3_cpu"};
	C6_unc = mat_data{dict_keys == "C6_cpu"};
	temp_unc = mat_data{dict_keys == "temp_pkg_cpu"};
	tpw =  mat_data{dict_keys == obj.TPW};
	tp0 =  mat_data{dict_keys == obj.TP0};

	res = [];
	
	% HERE I assume that cores in a socket are contiguous
	for sock=1:obj.num_socket
	
		idx1 = 1 + (sock-1)*cores_num;
		idx2 = (sock)*cores_num;
	
		res{sock}.temp = temp(:,idx1:idx2);
		%res{sock}.C0 = C0(:,idx1:idx2);
		res{sock}.mperf = mperf(:,idx1:idx2);
		res{sock}.freqd = freqd(:,idx1:idx2);
		res{sock}.clk = clk(:,idx1:idx2);
		res{sock}.clk_ref = clk_ref(:,idx1:idx2);
		res{sock}.top_down = top_down(:,idx1:idx2);
	
		%res{sock}.cstate = c(:,idx1:idx2,:);
		res{sock}.C0 = c(:,idx1:idx2,1);
		res{sock}.C1 = c(:,idx1:idx2,2);
		res{sock}.C3 = c(:,idx1:idx2,3);
		res{sock}.C6 = c(:,idx1:idx2,4);
		res{sock}.C7 = c(:,idx1:idx2,5);
		res{sock}.total_cstate = total_cstate(:,idx1:idx2);
	
		%%% CPU TODO CHANGE
		res{sock}.C0_unc = 1 - ( (C2_unc(:,sock)./ tsc(:,1)) + ...
			(C3_unc(:,sock)./ tsc(:,1)) + (C6_unc(:,sock)./ tsc(:,1)) );
		res{sock}.temp_unc = temp_unc(:,sock);
		res{sock}.tpw = tpw(:,sock);
		res{sock}.tp0 = tp0(:,sock);

	end
else
	res = mat_data;
end