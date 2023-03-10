large_interval = (0, 1000)
small_intervals = [(200, 300), (400, 500), (600, 700)]

# Sort the small intervals by their start points
small_intervals = sorted(small_intervals)

# Create a list of intervals covering the full range of the large interval
full_range = [(large_interval[0], large_interval[1])]
new_intervals = []
for interval in small_intervals:
    # Remove any small intervals that are outside the large interval
    if interval[0] >= large_interval[0] and interval[1] <= large_interval[1]:
        # Split the full range into two parts at the end point of each small interval
        new_intervals.extend([(large_interval[0], interval[0]), (interval[1], large_interval[1])])
    # If the small interval overlaps with the previous interval, merge them
    elif interval[0] < full_range[-1][1]:
        full_range[-1] = (full_range[-1][0], max(full_range[-1][1], interval[1]))
    # Otherwise, add the small interval to the new intervals list
    else:
        new_intervals.append((interval[0], interval[1]))
    # Update the full range with the new intervals
    full_range = merge_intervals(full_range, new_intervals)

# The result is the union of all the intervals in the full range list
result = [(interval[0], interval[1]) for interval in full_range]


print (result)