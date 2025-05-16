classdef THID
	
	properties
		Nmc = 4;
		Nmp = 2;
		%Nc = [];
		test_sampling_time = 0.5; %seconds

		num_cores = 56;
		num_socket = 2;

		model_order = 2;
		%num_prev_time = 20;
		num_add_eq = 10;

		metrics_core_names = ["C0", "freqd", "C6", "temp", "C1"];
		LEARNING_COMB_MODES = ["ALL_MEAN", "ALL_SINGLE", "ACTIVITY_SEPARATED"];

		% CONFIGURATIONS:
		TEST_TYPE = ["PRBS", "WAVE", "CLLAMA", "OLLAMA", "STREAM", "HPL"];
		TEST_NAMES = ["NOV2_prbs_random", "NOV2_waves_random", "NOV2_cllama", "NOV2_ollama", "NOV2_stream", "NOV2_hpl"];
		TEST_CONF = ["4", "8", "16", "24", "32", "48", "56", "112"];
		datapy_format = ["pkl", "csv"];
		datapy_type = ["list of dict", "dict", "list", "numpy_arr"];
		
		DICT_KEYS = ["temp_core", "clk_curr_core", "clk_ref_core", "aperf_core", "mperf_core", ...
			"tsc_core", "ThermStat_core", "C1_core", "C3_core", "C6_core", "C7_core", ...
			"fp_arith_inst_retired.128b_packed_double_core", "fp_arith_inst_retired.128b_packed_single_core", ...
			"fp_arith_inst_retired.256b_packed_double_core", "fp_arith_inst_retired.256b_packed_single_core", ...
			"fp_arith_inst_retired.512b_packed_double_core", "fp_arith_inst_retired.512b_packed_single_core", ...
			"freq_aperf_mperf_core", "busy_core", "topdown_retiring_slots_core", "topdown_bad_spec_slots_core", ...
			"topdown_fe_bound_slots_core", "topdown_be_bound_slots_core", "topdown_metrics_core", "topdown_slots_core", ...
			"temp_pkg_cpu", "erg_pkg_cpu", "erg_pp0_cpu", "PLR_cpu", "powCTL_cpu", "C2_cpu", "C3_cpu", ...
			"C6_cpu", "erg_psys_cpu", "power_from_erg_pkg_cpu", "power_from_erg_pp0_cpu"];
		
		TEMP_CORE = "temp_core";
		APERF_CORE = "aperf_core";
		MPERF_CORE = "mperf_core";
		FREQ = "freq_aperf_mperf_core";
		TSC = "tsc_core";
		TOP_DOWN = "topdown_be_bound_slots_core";
		CLK_CORE = "clk_curr_core";
		CLK_REF_CORE = "clk_ref_core";
		C_CORE = ["mperf_core", "C1_core", "C3_core", "C6_core", "C7_core" ];
		TPW = "power_from_erg_pkg_cpu";
		TP0 = "power_from_erg_pp0_cpu";
	end

	properties(SetAccess=protected, GetAccess=public)
		Pidle = [];
		Tidle = []
	end

	methods
		function obj = THID(number_socket, number_cores)
			
			obj.num_socket = number_socket;
			obj.num_cores = number_cores;

			obj.Pidle = zeros(number_socket,1);
			obj.Tidle = zeros(number_socket,1);
		end

		[M, control_M] = GetNProcess(obj, TEST,TEST_CONF, ROUND, dfmat_idx, dtype_idx, cco, preproc)
		res = parser_script(obj, ismetric, test_idx, tconf_idx, tround, dfmat_idx, dtype_idx, path_to_pm_file)
		[mc, Nc] = accMetrics(obj, M, learning_comb)
		[] = metricsAnalyze(obj, M, show, cores_idx)
		function [mcore, mcpu] = compileMetrics(obj,Mcpu, Tconst)
		
			dim = size(Mcpu.temp);
		
			mcore = nan(dim(1), obj.Nmc, dim(2));
			mcpu = nan(dim(1), obj.Nmp);
			
			mcore(:,1,:) = Mcpu.temp - Tconst(1);
			mcore(:,2,:) = Mcpu.freqd .* Mcpu.C0 * 10^6;
			mcore(:,3,:) = Mcpu.C1;
			%mcore(:,4,:) = M.C3;
			mcore(:,4,:) = Mcpu.C6;
			%mcore(:,4,:) = M.C7;
			%mcore(:,5,:) = Mcpu.top_down*10e-9;
			
			mcpu(:,1) = Mcpu.temp_unc' - Tconst(2);
			mcpu(:,2) = Mcpu.C0_unc';
		end
		function DesignMatrix = buildModelMatrix(obj, mc, munc, N, Nc, offset)
		
			% Reshape A for matrix construction
			A_reshaped = zeros(N, obj.Nmc * Nc);
			for i = 1:N
				A_t = squeeze(mc(i,:, :)); % Extract A(:, :, t)
				A_reshaped(i, :) = A_t(:)'; % Flatten and store
			end
			% Construct the system of equations
			% POW = [A_reshaped, B_expanded, kron(eye(t), ones(obj.Nmc, 1))] * Theta
			if offset
				offmat = ones(N,1)*obj.num_cores;
			else
				offmat = [];
			end
			DesignMatrix = [offmat, A_reshaped, munc];
		end
		function [Pi, Ti, obj] = computePidle(obj, M, c0_th_lvl, show)

			if ((nargin < 4) || isempty(show) )
				show = 0;
			end
			% Constant way
			Pi = nan(obj.num_socket,1);
			Ti = nan(obj.num_socket,1);
			% All but c0
			Pidle_stor = [];
			for cpu=1:obj.num_socket
				% start with all low
				idx = sum(M{cpu}.C0 < c0_th_lvl,2) == obj.num_cores;
				Pidle_stor{cpu}(:,1) = M{cpu}.tpw(idx);
				Pidle_stor{cpu}(:,2) = prctile(M{cpu}.temp(idx,:),80,2);
				%Pidle_stor{cpu}(:,2) = median(M{cpu}.temp(idx,:),2);
				%Pidle_stor{cpu}(:,2) = min(M{cpu}.temp(idx,:),[],2);
				%Pidle_stor{cpu}(:,2) = max(M{cpu}.temp(idx,:),[],2);

				if show
					figure();
					subplot(2,1,1);
					ax1 = plot(Pidle_stor{cpu}(:,1)); grid on;
					ylabel("Power [W]");
					subplot(2,1,2);
					ax2 = plot(Pidle_stor{cpu}(:,2), 'r'); grid on;
					ylabel("Power [W]");
					ax = findall(gcf, 'type', 'axes'); % Find all axes handles
					linkaxes(ax, 'x');
					pause(0.5);
				end
			end
			for cpu=1:obj.num_socket
				P = Pidle_stor{cpu}(:,1);
				T = Pidle_stor{cpu}(:,2);
				% Process data to just take a constant value:
				idx = (T < prctile(T,80));
				Ppr = P(idx);
				Tpr = T(idx);
				idx = (Ppr < prctile(Ppr,80));
				Ppr = Ppr(idx);
				Tpr = Tpr(idx);
				
				Pi(cpu) = mean(Ppr);
				Ti(cpu) = mean(Tpr);
			end

			obj.Pidle = Pi;
			obj.Tidle = Ti;
		end
		function [A,B,C] = getSSmat(obj, theta)
			A = nan(obj.model_order);
			A(:,1) = -theta(1:obj.model_order);
			A(1:end-1, 2:end) = eye(obj.model_order-1);
			A(end, 2:end) = 0;

			dim2 = (length(theta)-obj.model_order) / obj.model_order;
		
			Bdim = [obj.model_order, dim2]; %obj.num_cores+1];
			B = nan(Bdim);
			%{
			for mn = 1:obj.model_order
				idx1 = obj.model_order + (mn-1)*(Bdim(2)) + 1;
				idx2 = obj.model_order + (mn)*(Bdim(2));
				B(mn,:) = theta(idx1:idx2);
			end
			%}
			for mn = 1:Bdim(2)
				idx1 = obj.model_order + (mn-1)*(Bdim(1)) + 1;
				idx2 = obj.model_order + (mn)*(Bdim(1));
				B(:,mn) = theta(idx1:idx2);
			end
		
			C = zeros(1, obj.model_order);
			C(1) = 1;
		end
	
	end

	methods(Static)
		%[p,n]=numSubplots(n)
		function is_usable = isUsableNumber(x)
    		is_usable = ~isempty(x) && isnumeric(x) && isfinite(x);
			is_usable = is_usable & (x>=0) & (x<10^11);
		end		
		function or = fixData(M, fields, num_socket)
			% 1) Find the minimum among length
			for sk = 1:num_socket %manage more than 1 socket
				field = fields{1};
				minsz = length(M{sk}.(field));
				for i=1:numel(fields)
					field = fields{i};
					ttsz = length(M{sk}.(field));
					if (ttsz < minsz)
						minsz = ttsz;
					end
				end
				% 2) allign, then 3) remove firsts and lasts row:
				for i=1:numel(fields)
					field = fields{i};
					M{sk}.(field) = M{sk}.(field)(1:minsz,:);
					M{sk}.(field) = M{sk}.(field)(3:end-2,:);
				end
			end			
			or = M;
		end
		function U = computePower(mc, munc, pw_c, pw_unc)
			
			mc = permute(mc, [3 2 1]);

			%muly = pagemtimes(mc, pw_c);
			unc = munc*pw_unc(:);
			
			%U = squeeze(arrayfun(@(k) diag(muly(:,:,k)), 1:size(muly,3), 'UniformOutput', false));
			%U = cell2mat(U).'; % Transpose to get MxP

			for i=1:size(mc,3)
				U(i,:) = sum(mc(:,:,i).*pw_c',2)';
			end

			U = [U, unc];
		end
		function [Yn, Un] = filterCores(Y, U, core, Nrows, Ncols)

			% Here assuming cores are in an MxN rectangle with horizontal
			%	numeration
			crow = fix((core - 1) / Ncols) + 1;
			ccol = mod((core-1), Ncols) + 1;

			nrow = crow - 1;
			if nrow > 0
				North = (nrow-1)*Ncols + ccol;
			else
				North = -1;
			end
			nrow = crow + 1;
			if nrow > Nrows
				South = -1;
			else
				South = (nrow-1)*Ncols + ccol;
			end
			ncol = (ccol+1);
			if ncol > Ncols
				East = -1;
			else
				East = (crow-1)*Ncols + ncol;
			end
			ncol = (ccol-1);
			if ncol > 0
				West = (crow-1)*Ncols + ncol;
			else
				West = -1;
			end

			core_idx = [North, East, South, West];
			core_idx = core_idx(core_idx>0);

			Yn = Y(:,core);
			ll = size(U,2);
			Un = U(:, [core_idx ll]);
        end
        function [Yn, Un] = filterCoresRnd(Y, U, core, Nrows, Ncols, Nneigh)
		    % 	% Here assuming cores are in an MxN rectangle with horizontal
            min_id = 1;
            max_id = Nrows*Ncols; %first id is 1
            possibleValues = setdiff(min_id:max_id, core);
            core_idx = datasample(possibleValues, Nneigh, 'Replace', false);
            Yn = Y(:,core);
            ll = size(U,2);
            Un =  U(:, [core_idx ll]);
        end
    end
end