% Basic ORGBRAND (soft detection) 

% As introduced in
% K. R. Duffy, “Ordered reliability bits guessing random additive noise 
% decoding," in IEEE ICASSP, 2021, pp. 8268–8272 
% and implemented using the Landslide algorithm introduced in 
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
% random additive noise decoding,” arXiv:2202.13951, 2022.

% Inputs:
%   code_class      - CRC, else assumed defined by matrix
%   n               - code length
%   H               - Parity check matrix or CRC function
%   max_query       - Maximum number of code-book queries to abandonment
%   y_soft          - Channel soft information
%
% Outputs:
%   y_decoded       - Decoded ML codeword
%   err_vec         - Putative noise
%   n_guesses       - Number of guesses performed
%   abandoned       - 1 if abandoned, 0 if found a codeword

function [y_decoded,err_vec,n_guesses,abandoned] = bin_ORBGRAND(code_class,n,H,max_query,y_soft)

    % Hard demodulate
    y_demod = (y_soft<0);

    n_guesses = 1;
    err_vec = zeros(1,n);

    % Logistic starting weight
    W=0;

    %Abandon defualt values
    y_decoded = -1*ones(size(y_demod));

    % First query is demodulated string

    if isequal(code_class,'CRC')
        % MATLAB expects codeswords to be columns - H here is a CRC
        % check function
        t = mod(y_demod,2);
        [~,err] = H(mod(t',2));
        % If no error
        if ~err
            y_decoded = t;
            abandoned = 0;
            return;
        end
    % If not a CRC, assume you have the H matrix
    else 
        Hy = mod(H*y_demod',2);
        if Hy==zeros(size(Hy))
            y_decoded = y_demod;
            abandoned = 0;
            return;
        end
    end
    
    % Bit reliability
    reliability = abs(y_soft);
    [~, ind_order] = sort(reliability,'ascend');

    % Inverse sort order
    inv_perm = zeros(1,length(ind_order));
    for ii=1:length(ind_order)
           inv_perm(ind_order(ii))=ii;    
    end

    % If CRC
    if isequal(code_class,'CRC')
        while n_guesses<max_query
            % Increment logistic weight
            W=W+1;
            % Reset the Hamming weight.
            w = max(0,floor( ((1+2*n)-sqrt((1+2*n)^2-8*W))/2));
            while w<floor((sqrt(1+8*W)-1)/2)
                % Increment Hamming weight
                w=w+1;
                % Make error vectors
                noise_locations = landslide(W,w,n);
                % For each error vector
                for jj=1:size(noise_locations,1)
                    n_guesses = n_guesses +1;
                    err_vec =zeros(1,n);
                    err_vec(noise_locations(jj,:))=1;
                    err_vec = err_vec(inv_perm);
                    t = mod(y_demod-err_vec,2);
                    % MATLAB expects codeswords to be columns. H here is 
                    % a CRC check function
                    [~,err] = H(mod(t',2));
                    % If no error
                    if ~err
                        y_decoded = t;
                        abandoned = 0;
                        return;
                    end
                end
            end   
        end
    % If H a matrix we're dealing with a linear code
    elseif ismatrix(H)
        % This is the syndrome
        Hy = mod(H*y_demod',2);
        % This is the H columns reordered to put in ML order
        test_H = H(:,ind_order);
         while n_guesses<max_query
            % Increment logistic weight
            W = W+1;
            % Reset Hamming weight
            w = max(0,floor( ((1+2*n)-sqrt((1+2*n)^2-8*W))/2));
            % Increment Hamming weight
            while w<floor((sqrt(1+8*W)-1)/2)
                w=w+1;
                % Make error vectors
                % Internally converts W and n to W' and n'.
                noise_locations = landslide(W,w,n); 
                % For each error vector
                for jj=1:size(noise_locations,1)
                   n_guesses = n_guesses +1;
                    err_vec =zeros(1,n);
                    err_vec(noise_locations(jj,:))=1;
                    if (Hy == mod(test_H*err_vec',2))
                        err_vec = err_vec(inv_perm);
                        y_decoded = mod(y_demod-err_vec,2);
                        abandoned = 0;
                        return;
                    end
                end
            end
         end
    end

    % If we max out on queries.
    abandoned = 1;
    err_vec = zeros(size(y_demod));

end

