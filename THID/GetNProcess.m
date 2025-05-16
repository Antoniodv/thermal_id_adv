function [M, control_M] = GetNProcess(obj, TEST,TEST_CONF, ROUND, dfmat_idx, dtype_idx, cco, preproc)
	M = [];
	fields = [];
	%control group
	control_M = [];

	for test_idx = TEST
		for tci=TEST_CONF
			for tri=ROUND
				res = obj.parser_script(1, test_idx, tci, tri, dfmat_idx, dtype_idx);
				if ~isempty(res)
					if isempty(M) % first populating M and fields
						fields = fieldnames(res{1});
						M = res;
						% pre-processing
						%	- align all the fields, since they have different size..... 
						%	- remove firsts and lasts rows since data are bad
						M = obj.fixData(M, fields, obj.num_socket);
					else % not empty M
						if (tci == cco(1)) && (tri == cco(2)) %setup control group 
							controlM = res;
							%removing firsts and lasts rows since data are bad
							controlM = obj.fixData(controlM, fields, obj.num_socket);
						else % not control group
							res = obj.fixData(res, fields, obj.num_socket);
							for sk = 1:obj.num_socket
								for i=1:numel(fields)
									field = fields{i};
									M{sk}.(field) = [M{sk}.(field); res{sk}.(field)];
								end
							end
						end
					end
				end
			end
		end
	end

	if preproc
		nummiss = 0;
		missTable = table('Size', [0, 5], 'VariableTypes', {'int32','int32','int32', 'string','double'}, ...
                  		'VariableNames', {'row', 'column', 'cpu', 'field', 'value'});		
		for sk=1:obj.num_socket
			for i=1:numel(fields)
				field = fields{i};
				lenfin = size(M{sk}.(field),1);
				r = 1;
				while (r<=lenfin)
					for c=1:size(M{sk}.(field),2)
						if ~obj.isUsableNumber( M{sk}.(field)(r,c) )
							nummiss = nummiss + 1;
							missTable(nummiss,:) = {r, c, sk, field, M{sk}.(field)(r,c)};
							savefield = field;
							for i2=1:numel(fields)
								field = fields{i2};
								M{sk}.(field)(r,:) = [];
							end
							field = savefield;
							lenfin = lenfin-1;
							r = r-1;
							break;
						end
					end
					r = r+1;
				end
			end
		end
		disp(['Number of problematic values: ', num2str(nummiss)]);
		disp(missTable);
	end
end