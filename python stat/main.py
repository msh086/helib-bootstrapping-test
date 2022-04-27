import re
import numpy as np


# return (basic info, [sum(mean), sum(var), sum(mag) for every level], sample_count)
def collect_data(fname, log2q):
    with open(fname) as f:
        header = f.readline()[:-1]  # remove the '\n' char
        param_idx, total_levels, boot_levels, square_between_boots = [int(_) for _ in header.split()]
        lines_per_sample = total_levels + 1 + boot_levels + square_between_boots * (boot_levels - 1)
        # sum(x)/dim, sum(x^2)/dim, sum(mag) for every level
        # TODO: put mag into a list to get xx% confidence interval for mag?
        stat = np.zeros([lines_per_sample, 3], float)
        sample_stat = np.zeros_like(stat)
        # modulus size for every level
        mod_stat = [0. for _ in range(lines_per_sample)]
        sample_count = 0
        # read a whole sample each time
        while True:
            line = f.readline()
            match = re.match(r"SAMPLE \d+", line)
            if match is None:
                break
            for i in range(lines_per_sample):
                line = f.readline()
                if len(line) == 0:
                    return (param_idx, total_levels, boot_levels, square_between_boots, mod_stat), stat, sample_count
                line = line[:-1]
                parts = line.split(', ')
                if len(parts) != 6:
                    print("Error: format")
                    exit(1)
                cap, noise, mod, mean, var, mag = [float(_) for _ in parts]
                if mod_stat[i] == 0:
                    mod_stat[i] = mod
                elif mod_stat[i] != mod:  # debug
                    print("Warning: {}-th sample, {}-th line, old mod: {}, new mod: {}"
                          .format(sample_count, i, mod_stat[i], mod))
                mean = abs(mean)  # NOTE: always positive...
                # normalize
                dist = log2q - mod
                mean += dist
                var += dist * 2
                mag += dist
                sample_stat[i][0] += 2.0 ** mean
                sample_stat[i][1] += 2.0 ** var + (2.0 ** mean) ** 2
                sample_stat[i][2] += 2.0 ** mag
            # a whole sample is read, update stat
            sample_count += 1
            stat += sample_stat
        return (param_idx, total_levels, boot_levels, square_between_boots, mod_stat), stat, sample_count


def merge_stats(stat_list):
    nfiles = len(stat_list)
    if nfiles == 0:
        print("Error: empty file list")
        exit(1)
    config, stat, sample_count = stat_list[0]
    for _config, _stat, _sample_count in stat_list[1:]:
        if _config != config:
            print("Error: config mismatch")
            print("Expecting: ", config)
            print("Got: ", _config)
            exit(1)
        sample_count += _sample_count
        stat += _stat
    stat[:, 0] /= sample_count  # E(x)
    stat[:, 1] /= sample_count  # E(x^2)
    stat[:, 2] /= sample_count  # E(mag)
    stat = np.log2(np.hstack([stat, (stat[:, 1] - stat[:, 0]**2)[:, np.newaxis]]))
    return config, stat, sample_count


data_files = ["output_{}".format(i) for i in range(8)]

if __name__ == '__main__':
    stat_list = []
    log2q = 551.028904
    for fname in data_files:
        stat_list.append(collect_data(fname, log2q))
    config, stat, sample_count = merge_stats(stat_list)
    print("number of samples: {}".format(sample_count))
    print("log2q = {}".format(log2q))
    print("mean var mag ")
    slicer = [0, 3, 2]
    print(stat[:, slicer])
