#! /usr/bin/env python

reads = [[0, 76], [7, 83], [22, 98], [39, 115], [44, 120], [71, 147], [76, 152], [79, 155], [83, 159], [86, 162], [88, 164], [90, 166], [91, 167], [102, 178], [102, 178], [117, 193], [122, 198], [126, 202], [126, 202], [126, 202], [130, 206], [131, 207], [132, 208], [132, 208], [135, 211], [136, 212], [136, 212], [169, 245], [172, 248], [185, 261], [192, 268], [202, 278], [209, 285], [221, 297], [224, 300], [235, 311], [264, 340], [279, 355], [286, 362], [297, 373], [301, 377], [303, 379], [307, 383], [310, 386], [313, 389], [320, 396], [322, 398], [325, 401], [326, 402], [333, 409], [340, 416], [350, 426], [360, 436], [399, 475], [410, 486], [411, 487], [412, 488], [421, 497], [434, 510], [440, 516], [458, 534], [538, 614], [541, 617], [541, 617], [547, 623], [549, 625], [553, 629], [584, 660]]

pairedLens = {384: 1, 386: 1, 261: 1, 326: 2, 328: 1, 172: 1, 458: 1, 466: 1, 219: 1, 412: 1, 350: 1, 288: 1, 99: 1, 229: 1, 166: 1, 295: 1, 300: 2, 301: 1, 432: 1, 245: 1, 251: 1, 202: 1, 319: 1}

unique_reads = [reads[0]]
read_counts = [1]
for r in reads[1:]:
    if r == unique_reads[-1]:
        read_counts[-1] += 1
    else:
        unique_reads.append(r)
        read_counts.append(1)

matches = dict()
for i in range(len(unique_reads)-1):
    for j in range(i+1, len(unique_reads)):
        l = unique_reads[j][1] - unique_reads[i][0]
        if l in pairedLens:
            if l in matches:
                matches[l].append((i,j))
            else:
                matches[l] = [(i,j)]

matches = sorted(matches.items(), key=lambda x : len(x[1]))
for m in matches:
    print(str(m[0]) + '\t' + str(m[1]))

num_matches = len(matches)
lenCounts = [0] * num_matches
for i in range(num_matches):
    lenCounts[i] = pairedLens[matches[i][0]]

paths = set()
total_paths = 0

def solve_recursive(matches_id, counts, depth=0, path=[]):
    global total_paths
    if matches_id >= num_matches:
        total_paths += 1
        #paths.add(str(path))
        #unpaired = []
        #for i in range(len(unique_reads)):
        #    for _ in range(counts[i]):
        #        unpaired.append(unique_reads[i])
        #a = [[[], unpaired]]
        #return a
        return

    lenCounts[matches_id] -= 1
    if lenCounts[matches_id] == 0:
        new_matches_id = matches_id+1
    else:
        new_matches_id = matches_id

    #results = []

    for p in matches[matches_id][1]:
        if counts[p[0]] > 0 and counts[p[1]] > 0:
            counts[p[0]] -= 1
            counts[p[1]] -= 1

            solve_recursive(new_matches_id, counts, depth+1, path+[p])

            counts[p[0]] += 1
            counts[p[1]] += 1

            #if len(results) > 0:
            #    for r in res:
            #        r[0].append([unique_reads[p[0]], unique_reads[p[1]]])
            #    results += res

    lenCounts[matches_id] += 1

    #return results

solve_recursive(0, read_counts)
#print(len(paths))
print(total_paths)