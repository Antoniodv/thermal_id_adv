function [A,B,C] = getSSmat(theta, model_order)
			A = nan(model_order);
			A(:,1) = -theta(1:model_order);
			A(1:end-1, 2:end) = eye(model_order-1);
			A(end, 2:end) = 0;

			dim2 = (length(theta)-model_order) / model_order;
		
			Bdim = [model_order, dim2]; %obj.num_cores+1];
			B = nan(Bdim);
			%{
			for mn = 1:model_order
				idx1 = model_order + (mn-1)*(Bdim(2)) + 1;
				idx2 = model_order + (mn)*(Bdim(2));
				B(mn,:) = theta(idx1:idx2);
			end
			%}
			for mn = 1:Bdim(2)
				idx1 = model_order + (mn-1)*(Bdim(1)) + 1;
				idx2 = model_order + (mn)*(Bdim(1));
				B(:,mn) = theta(idx1:idx2);
			end
		
			C = zeros(1, model_order);
			C(1) = 1;
end