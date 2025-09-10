function subplotlabel(char, location, fontsize)
    % Set default location if not provided
    if nargin < 2 || isempty(location)
        location = 'northeast';  % Default location
    end
    
    % Set default font size if not provided
    if nargin < 3 || isempty(fontsize)
        fontsize = 12;  % Default font size
    end

    % Create the label with bold formatting
    lbl = sprintf('\\bf{%s}', char);

    % Determine the text position and padding based on the location
    switch location
        case 'northwest'
            xPos = 0; yPos = 1;
            hAlign = 'left'; vAlign = 'top';
            lblText = [' ' lbl];  % Add space to the left of the label
        case 'northeast'
            xPos = 1; yPos = 1;
            hAlign = 'right'; vAlign = 'top';
            lblText = [lbl ' '];  % Add space to the right of the label
        case 'southwest'
            xPos = 0; yPos = 0;
            hAlign = 'left'; vAlign = 'bottom';
            lblText = [' ' lbl];  % Add space to the left of the label
        case 'southeast'
            xPos = 1; yPos = 0;
            hAlign = 'right'; vAlign = 'bottom';
            lblText = [lbl ' '];  % Add space to the right of the label
        otherwise
            error('Invalid location. Choose from northwest, northeast, southwest, or southeast.');
    end

    % Add the text to the subplot
    text(xPos, yPos, lblText, 'Units', 'normalized', 'VerticalAlignment', vAlign, 'HorizontalAlignment', hAlign, 'FontSize', fontsize);
end
