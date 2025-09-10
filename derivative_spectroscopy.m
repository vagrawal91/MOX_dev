function [A, dA, d2A, valid_lambda] = derivative_spectroscopy(A0, lambda, filter_scale, order)
    % Function to calculate smoothed spectra and their first and second
    % derivatives using the Savitzky-Golay filter.
    %
    % Inputs:
    %   A0           - Matrix of original spectra (rows = individual spectra, columns = wavelengths).
    %   lambda       - Wavelength vector (assumed equidistant).
    %   filter_scale - Filter wavelength scale for smoothing.
    %   order        - Polynomial order of the Savitzky-Golay filter.
    %
    % Outputs:
    %   A            - Smoothed spectra matrix for valid indices.
    %   dA           - First derivative matrix for valid indices.
    %   d2A          - Second derivative matrix for valid indices.
    %   valid_lambda - Wavelength vector for valid indices.

    % Check input dimensions
    [num_samples, num_wavelengths] = size(A0);
    if length(lambda) ~= num_wavelengths
        error('The number of columns in A0 must match the length of lambda.');
    end

    % Calculate the sampling interval in the wavelength domain
    delta_lambda = mean(diff(lambda));

    % Compute window size based on filter scale
    filter_width = round(filter_scale / delta_lambda); % Convert filter scale to points
    if mod(filter_width, 2) == 0
        filter_width = filter_width + 1; % Ensure window size is odd
    end

    % Validate polynomial order
    if order >= filter_width
        error('Polynomial order (%d) must be less than the filter width (%d).', order, filter_width);
    end
    if order < 2
        error('Polynomial order (%d) must be at least 2 for the second derivative to exist.', order);
    end

    % Preallocate outputs
    A = [];
    dA = [];
    d2A = [];

    % Calculate Savitzky-Golay filter coefficients
    [~, g] = sgolay(order, filter_width);

    % Loop through each spectrum
    for i = 1:num_samples
        % Apply smoothing with Savitzky-Golay filter
        Smoothed = conv(A0(i, :), g(:, 1), 'same');

        % Compute the first derivative
        dA_row = conv(A0(i, :), -1/delta_lambda * g(:, 2), 'same');

        % Compute the second derivative
        d2A_row = conv(A0(i, :), 2/(delta_lambda^2) * g(:, 3), 'same');

        % Determine valid indices
        half_width = (filter_width - 1) / 2;
        valid_indices = false(1, num_wavelengths);
        valid_indices(half_width+1:end-half_width) = true;

        % Restrict results to valid indices
        A = [A; Smoothed(valid_indices)];
        dA = [dA; dA_row(valid_indices)];
        d2A = [d2A; d2A_row(valid_indices)];
        valid_lambda = lambda(valid_indices); % Update for each iteration
    end
end
