import sys
import os
import argparse
import fire


def make_dep_dict(fai, depth):
    a = {}
    with open(fai, 'r') as f:
        for line in f:
            line = line.rstrip().split()
            chr_id = line[0]
            chr_len = int(line[1])
            tmp = [0 for x in range(chr_len)]
            with open(depth, 'r') as g:
                for info in g:
                    info = info.rstrip().split()
                    if info[0] == chr_id:
                        pos = int(info[1]) - 1
                        tmp[pos] = float(info[2])
            a[chr_id] = tmp

    return a


def fire_dict(input_fire, fai, bin_len, depth):
    depth_dict = make_dep_dict(fai, depth)
    bin_len = int(bin_len)
    a = {}
    with open(input_fire, 'r') as f:
        for line in f:
            line = line.rstrip().split()
            range_bin = str(line[0]) + "_" + str(line[1]) + "_" + str(line[2])
            a[range_bin] = line[3]

    with open(fai, 'r') as g:
        for info in g:
            info = info.rstrip().split()
            chr_id = info[0]
            chr_len = int(info[1])
            z = open(str(chr_id) + "-" + str(bin_len) + ".fire", 'w')
            out = open(str(chr_id) + "-" + str(bin_len) + ".rna.depth", 'w')
            for x in range(chr_len / int(bin_len) + 1):
                start = x * bin_len
                end = (x + 1) * bin_len
                if end > chr_len:
                    end = chr_len
                tmp_range_bin = str(chr_id) + "_" + str(start) + "_" + str(end)

                if tmp_range_bin in a.keys():
                    value = a[tmp_range_bin]
                else:
                    value = 0.000001

                if end < chr_len:
                    tmp_depth = sum(depth_dict[chr_id][start:end])
                else:
                    tmp_depth = sum(depth_dict[chr_id][start:])
                mid = start + bin_len / 2
                z.write(str(mid) + '\t' + str(value) + '\n')
                out.write(str(mid) + '\t' + str(tmp_depth) + '\n')

            z.close()
            out.close()


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fai", type=str, help="the chr length file")
    parser.add_argument("-r", "--res", type=int, help="the resolution of the fire file")
    parser.add_argument("-d", "--depth", type=str, help="the RNA depth file")
    parser.add_argument("-m", "--fire", type=str, help="the fire value file")
    command = sys.argv[1:]
    args = parser.parse_args(command)
    return command, args


def main():
    command, args = get_args()
    fai = args.fai
    res = args.res
    depth = args.depth
    fire = args.fire
    fire_dict(fire, fai, res, depth)


if __name__ == "__main__":
    main()
