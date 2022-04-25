import re


def collect_data(fname):
    with open(fname) as f:
        header = f.readline()[:-1]  # remove the '\n' char
        param_idx, total_levels, boot_levels, square_between_boots = [int(_) for _ in header.split()]
        lines_per_sample = total_levels + 1 + boot_levels + square_between_boots * (boot_levels - 1)
        # sum(mean), sum(var), sum(mag) for every level
        stat = [[0 for _ in range(3)] for _ in range(lines_per_sample)]
        # modulus size for every level
        mod_stat = [0 for _ in range(lines_per_sample)]
        lines_remaining = 0
        sample_count = 0
        while True:
            line = f.readline()
            if len(line) == 0:
                return (param_idx, total_levels, boot_levels, square_between_boots, sample_count), stat
            line = line[:-1]
            if lines_remaining == 0:
                match = re.match(r"SAMPLE \d+", line)
                if match is None:
                    exit(1)
                lines_remaining = lines_per_sample
            else:
                parts = line.split(', ')
                if len(parts) != 6:
                    exit(1)


data_files = []


if __name__ == '__main__':
    pass