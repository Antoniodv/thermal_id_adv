function theta_vector_remapped = remap_theta_vector(theta_vector, num_cores, active_cores_indexes)

        %exogen
        exogen_tmp = theta_vector(3:end);

        % do we have enough exo coeff to cover all active cores?
        if length(exogen_tmp) < 2*length(active_cores_indexes)
            error("not enough exo coefficients to cover all active cores")
        end

        % do we have too many exo coeff to cover all active cores?
        if length(exogen_tmp) > 2*length(active_cores_indexes)
            error("too many exo coefficients to cover all active cores")
        end


        
        % any core index over the maximum value allowed by the total
        % number of cores in the whole processor?
        if any(active_cores_indexes > num_cores)
            error('at least one core index is over the maximum allowed from num_cores');
        end
        
        if any(active_cores_indexes < 1)
            error('at least one core index is lower than 1');
        end

        theta_vector_remapped_exo = zeros(1, 2*num_cores);
        for j = 1:length(active_cores_indexes)
         core_index = active_cores_indexes(j);
         global_index = (2*core_index)-1;
         
         exo_index = (j-1)*2 + 1;
         
         theta_vector_remapped_exo(global_index) = exogen_tmp(exo_index);
         theta_vector_remapped_exo(global_index+1) = exogen_tmp(exo_index+1);
        end
        theta_vector_remapped = [ theta_vector(1:2), theta_vector_remapped_exo ];

end