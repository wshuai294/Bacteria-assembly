

window = 1000
repeat_num = 0
count_list = []

for line in open("/home/wangshuai/assembly_result/real_data_w2/DRR198803.ref.origin.kmer.count.tab"):
    array = line.strip().split()
    # locus = int(array[0])
    count = int(array[1])
    if count > 1:
        count_list.append(1)
    else:
        count_list.append(0)


start, end = 0, 0
continuous = False


uniq_interval = []
for i in range(len(count_list) - window + 1):
    repeat_num = sum(count_list[i:i+window])
    repeat_ratio = repeat_num/window
    

    if repeat_ratio < 0.1:

        if continuous == False:
            start = i
            end = i + window
        else:
            end = i + window
        continuous = True
        # print (i, repeat_ratio)
    else:
        if continuous == True:
            # print (start, end)
            if len(uniq_interval) > 0 and start <= uniq_interval[-1][1]:
                uniq_interval[-1][1] = end
            else:
                uniq_interval.append([start, end])
        continuous = False

if continuous == True:
    if len(uniq_interval) > 0 and start <= uniq_interval[-1][1]:
        uniq_interval[-1][1] = end
    else:
        uniq_interval.append([start, end])

uniq_len = 0
for interval in uniq_interval:
    interval_len = (interval[1] - interval[0])
    if interval_len < 10000:
        continue
    print (interval)
    uniq_len += interval_len

print (uniq_len)
print (count_list[2510352:2512217])

