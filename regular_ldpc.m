function H = regular_ldpc(message_length, codeword_length, num_of_ones)
    % Initialize the parity check matrix
    H = zeros(codeword_length - message_length, codeword_length);

    % Fill the parity check matrix
    for i = 1:codeword_length
        % Get the indices for the "1"s in the column
        ones_indices = randperm(codeword_length - message_length, num_of_ones);
        
        % Set the "1"s in the column of the parity check matrix
        H(ones_indices, i) = 1;
    end

    % Adjust the matrix until it becomes a regular LDPC code
    while any(sum(H, 2) ~= num_of_ones)
        for i = 1:size(H, 1)
            while sum(H(i, :)) > num_of_ones
                % Find a "1" in the row
                one_index = find(H(i, :) == 1, 1);
                    
                % Find a "0" in the row
                zero_index = find(H(i, :) == 0, 1);
                
                % Swap the "1" and the "0"
                H(i, one_index) = 0;
                H(i, zero_index) = 1;
            end
        end
    end
end

