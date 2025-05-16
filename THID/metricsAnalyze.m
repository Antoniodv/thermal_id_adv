function [] = metricsAnalyze(obj, M, show, cores_idx)

	% power_from_erg_pp0_cpu_nan_avg is empty
	if show(1)
		for cpu=1:obj.num_socket
			Np = length(M{cpu}.tpw);
			t = 0:obj.test_sampling_time:(Np-1)*obj.test_sampling_time;
			figure();
			plot(t,M{cpu}.tpw); hold on; grid on;
			plot(t,M{cpu}.tp0);
			legend(["tpw", "tp0"]);
			ylabel("Power [W]");
			xlabel("Time [S]");
			pause(0.5);
		end
	end

	% Check C0 Uncore vs core
	if show(2)
		for cpu=1:obj.num_socket
			Np = length(M{cpu}.tpw);
			t = 0:obj.test_sampling_time:(Np-1)*obj.test_sampling_time;
			figure();
			plot(t,M{cpu}.C0_unc); hold on; grid on;
			plot(t,min(M{cpu}.C0,[],2), "Linewidth", 2);
			plot(t,max(M{cpu}.C0,[],2), "Linewidth", 2);
			legend(["C0 Uncore", "min(C0 cores)", "Max(C0 cores)"]);
			ylabel("[%]");
			xlabel("Time [S]");
		end

		[p,n] = numSubplots(length(cores_idx));
		for cpu=1:obj.num_socket
			Np = length(M{cpu}.tpw);
			t = 0:obj.test_sampling_time:(Np-1)*obj.test_sampling_time;
			figure(); title("sum of all c states = 1?");
			plotidx = 1;
			for core=cores_idx
				subplot(p(1),p(2),plotidx);
				plotidx = 1+plotidx;
				acc = M{cpu}.C0(:,core) + M{cpu}.C1(:,core) + ...
					M{cpu}.C3(:,core) + M{cpu}.C6(:,core) + M{cpu}.C7(:,core);
				plot(t,acc, "Linewidth", 1); hold on; grid on;
				title(int2str(core));
			end
			pause(0.5);
		end
	end

	% Check core C0 vs rest:
	if show(3)
		[p,n] = numSubplots(length(cores_idx));
		for cpu=1:obj.num_socket
			Np = length(M{cpu}.tpw);
			t = 0:obj.test_sampling_time:(Np-1)*obj.test_sampling_time;
			figure();
			plotidx = 1;
			for core=cores_idx
				subplot(p(1),p(2),plotidx);
				plotidx = 1+plotidx;
				plot(t,M{cpu}.C0(:,core), "Linewidth", 3); hold on; grid on;
				plot(t,M{cpu}.C1(:,core));
				plot(t,M{cpu}.C3(:,core));
				plot(t,M{cpu}.C6(:,core));
				plot(t,M{cpu}.C7(:,core), "Linewidth", 2);
				title(int2str(core));
			end
			legend(["C0", "C1", "C3", "C6", "C7"])
			ylabel("[%]");
			xlabel("Time [S]");
			pause(0.5);
		end
	end

	%Data dependency:
	if show(4)
		[test, Nc] = obj.accMetrics(M,"ALL_MEAN");
		for cpu=1:obj.num_socket
			test{cpu}.temp_unc = M{cpu}.temp_unc;
			test{cpu}.C0_unc = M{cpu}.C0_unc;
		end
		for cpu=1:obj.num_socket
			Ntot = length(test{cpu}.temp);			

			[metrics, metrics_cpu] = obj.compileMetrics(test{cpu}, [0 0]);
			
			DesignMatrix = obj.buildModelMatrix(metrics, metrics_cpu, Ntot, 1, 0);
			
			A = (M{cpu}.tpw - obj.Pidle(cpu));
			for i=1:(size(DesignMatrix,2)-1)
				figure();
				if i<obj.Nmc
					fact = obj.num_cores;
				else
					fact=1;
				end
				scatter(DesignMatrix(:,i)/fact,A); hold on; grid on;
				val_x = xlim;
				val_y = ylim;
				pp = [val_x; val_y];
				tpw = [DesignMatrix(:,i)/fact ones(length(DesignMatrix),1)] \ A;
				ll = val_x*tpw(1) + tpw(2);
				plot(val_x,ll);
				pause(0.5);
			end
		end
	end
	

	
end

