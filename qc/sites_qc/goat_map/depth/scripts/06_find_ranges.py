from pathlib import Path
import gzip


def is_within_thresholds(all_depth, low_depth, high_depth):
    if 'NA' in [all_depth, low_depth]:
        return False
    if high_depth == 'NA\n':
        return False

    all_depth = int(all_depth)
    if all_depth < 43 or all_depth > 376:
        return False

    low_depth = int(low_depth)
    if low_depth < 31 or low_depth > 278:
        return False

    high_depth = int(high_depth)
    if high_depth < 14 or high_depth > 106:
        return False

    return True


def find_ranges(input_file, output_file_in, output_file_out):
    in_file = open(output_file_in, 'w+')
    out_file = open(output_file_out, 'w+')
    is_in_range = False
    range_start_pos = 0
    range_end_pos = -1
    chromosome_name = '.'.join(input_file.stem.split('.')[:2])

    with gzip.open(input_file, 'rt') as f:
        first_line = True
        for line in f:
            # Ignore header.
            if first_line:
                first_line = False
                continue

            # Extract values.
            pos, all_depth, low_depth, high_depth = line.split('\t')
            pos = int(pos)-1

            if is_in_range:
                if is_within_thresholds(all_depth, low_depth, high_depth):
                    continue
                else:
                    range_end_pos = pos
                    is_in_range = False
                    in_file.write(f'{chromosome_name}\t{range_start_pos}\t{range_end_pos}\n')
            else:
                if is_within_thresholds(all_depth, low_depth, high_depth):
                    range_start_pos = pos
                    is_in_range = True
                    if range_end_pos == -1:
                        if range_start_pos == 0:
                            continue
                        out_file.write(f'{chromosome_name}\t{range_end_pos+1}\t{range_start_pos}\n')
                        continue
                    out_file.write(f'{chromosome_name}\t{range_end_pos}\t{range_start_pos}\n')

                else:
                    continue

        # End of file. Write final range.
        if is_in_range:
            in_file.write(f'{chromosome_name}\t{range_start_pos}\t{pos+1}\n')
        else:
            out_file.write(f'{chromosome_name}\t{range_end_pos}\t{pos+1}\n')

    in_file.close()


def main():
    work_dir = Path('/davidData/ffs/leopard/depth/output/combined_perChr/')
    for chromosome in work_dir.glob('*pos.gz'):
        if chromosome.parents[0].joinpath(f'{chromosome.stem}.range.in').is_file():
            continue
        find_ranges(chromosome,
                    chromosome.parents[0].joinpath(f'{chromosome.stem}.range.in'),
                    chromosome.parents[0].joinpath(f'{chromosome.stem}.range.out'))


if __name__ == '__main__':
    main()
