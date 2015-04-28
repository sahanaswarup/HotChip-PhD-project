function invDir = inverse_direction(direction)
switch direction
    case 'left'
        invDir = 'right';
    case 'right'
        invDir = 'left';
    case 'front' 
        invDir = 'back';
    case 'back'
        invDir = 'front';
    case 'up'
        invDir = 'down';
    case 'down'
        invDir = 'up';
        
end
