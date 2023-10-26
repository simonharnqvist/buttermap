def read_bed(file):
    with open(file, "r") as f:
        lines = f.readlines()

    regions = []
    for line in lines:
        chrom, start, end = line.split("\t")
        end = end.strip("\n")
        regions.append(f"{chrom}:{start}-{end}")

    return regions