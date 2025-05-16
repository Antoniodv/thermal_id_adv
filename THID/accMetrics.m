function [mc, Nc] = accMetrics(obj, M, learning_comb)
	mc = [];
	cthres = 0.6;
	fields = fieldnames(M{1});	
	
	switch learning_comb
		case obj.LEARNING_COMB_MODES(1) %"ALL_MEAN"
			Nc = 1;
			for cpu=1:obj.num_socket
				Ntot = length(M{cpu}.temp);
				for ee=1:Ntot
					for i=1:numel(fields)
						field = fields{i};
						if sum(field == obj.metrics_core_names)>0
							avgval = sum(M{cpu}.(field)(ee, :));
							if ~obj.isUsableNumber(avgval)
								avgval = 0;
								%end
							end
							mc{cpu}.(field)(ee,1) = avgval;
						end
					end
				end
			end
		case obj.LEARNING_COMB_MODES(2) %"ALL_SINGLE"
			Nc = obj.num_cores;
			mc = M;
		case obj.LEARNING_COMB_MODES(3) %"ACTIVITY_SEPARATED" --> merge active and inactive cores
			Nc = 3;
			
			for cpu=1:obj.num_socket
				Ntot = length(M{cpu}.temp);
				for ee=1:Ntot
				
					idxact = M{cpu}.C0(ee,:) > cthres;
					idxidle = M{cpu}.C1(ee,:) > cthres;
					idxothers = ~(idxact | idxidle);
				
					for i=1:numel(fields)
						field = fields{i};
						if sum(field == obj.metrics_core_names)>0
							avgval = sum(M{cpu}.(field)(ee, idxact));
							if ~obj.isUsableNumber(avgval)
								%{
								if sum(field == ["freqd", "temp"])>0
									if ee ~= 1 
										avgval = mc{cpu}.(field)(ee-1,1);
									else
										avgval = mean(M{cpu}.(field)(ee,:));
									end
								else
								%}
									avgval = 0;
								%end
							end
							mc{cpu}.(field)(ee,1) = avgval;
							avgval = sum(M{cpu}.(field)(ee, idxidle));
							if ~obj.isUsableNumber(avgval)
								%{
								if sum(field == ["freqd", "temp"])>0
									if ee ~= 1 
										avgval = mc{cpu}.(field)(ee-1,2);
									else
										avgval = mean(M{cpu}.(field)(ee,:));
									end
								else
								%}
									avgval = 0;
								%end
							end
							mc{cpu}.(field)(ee,2) = avgval;
							avgval = sum(M{cpu}.(field)(ee, idxothers));
							if ~obj.isUsableNumber(avgval)
							%{
								if sum(field == ["freqd", "temp"])>0
									if ee ~= 1 
										avgval = mc{cpu}.(field)(ee-1,3);
									else
										avgval = mean(M{cpu}.(field)(ee,:));
									end
								else
							%}
									avgval = 0;
								%end
							end
							mc{cpu}.(field)(ee,3) = avgval;
						end
					end
				end
			end %for cpu=1:obj.num_socket
		otherwise
			disp("Wrong cores combination type for learning!");
	end %switch

	%obj.Nc = Nc;
end