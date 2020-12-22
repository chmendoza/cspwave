function ind = add_buffer(ind, buffer_length, range)
% ind is an nx2 vector of n intervals in [range(1) range(2)]. ind(i,:) has 
% the start and end point of the i-th interval. buffer_length is the number 
% of samples to the right and to the left that the interval should be
% increased.

I = ind(:,1) - buffer_length >= range(1);
J = ind(:,2) + buffer_length <= range(2);
% If out of boundaries, use extreme points
ind(~I,1) = range(1); ind(~J,2) = range(2);
% Add buffer where there is no boundary violation
ind(I,1) = ind(I,1) - buffer_length;
ind(J,2) = ind(J,2) + buffer_length;

end