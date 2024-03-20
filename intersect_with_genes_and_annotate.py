import os, gzip, random
from collections import defaultdict


def extract_QC_passed_samples(QC_passed_samples):
    IDs_of_passed = {}
    with open(QC_passed_samples) as f:
        f.readline()
        for line in f:
            splitted_line = line.split("\t")
            IDs_of_passed[splitted_line[0]] = splitted_line[1]
    return(IDs_of_passed)

def extract_coords(gene_list_file, folders):
    prepared_coords = {}
    for i, gene_list in enumerate(gene_list_file):
        with open(gene_list) as f:
            genes_and_coords = defaultdict(list)
            for line in f:
                splitted_line = line.strip().split("\t")
                genes_and_coords["chr" + splitted_line[1]].append((splitted_line[5], splitted_line[6], splitted_line[1]))
        prepared_coords[folders[i]] = genes_and_coords
    return(prepared_coords)

def extract_cnvs(path_to_file):
    header = []
    cnvs = []
    with gzip.open(path_to_file, "r") as f:
        for line in f:
            line = line.decode("utf-8").encode("ascii")
            if line.startswith("##"):
                header.append(line.strip())
            elif line.startswith("#"):
                splitted_line = line.strip().split("\t")
                header.append(line.strip())
            else:
                splitted_line = line.strip().split("\t")
                if not splitted_line[0].startswith("chr"):
                    splitted_line[0] = "chr" + splitted_line[0]
                if True:
                    cnvs.append(splitted_line)
    return(header, cnvs)

def are_intersected(cnv, gene):
    x1 = int(cnv[1])
    x2 = int(cnv[2])
    y1 = int(gene[0])
    y2 = int(gene[1])
    return(x2 >= y1 and y2 >= x1)





def intersect_samples(list_of_gene_lists, results_folder, folders, path_to_output_folder, QC_passed_samples, calculate_frequencies, annotate_frequencies):
    random.seed(1)
    gene_lists = extract_coords(list_of_gene_lists, folders)
    number_of_SVs = {}
    names_of_samples = {}
    counter = 0.0
    QC_passed_samples = extract_QC_passed_samples(QC_passed_samples)


    if annotate_frequencies:
        frequencies = defaultdict(dict)
        for root, dirs, files in os.walk("/path/frequencies/"):
            for file in files:
                if file.endswith(".seg"):
                    with open("/path/frequencies/" + file) as f:
                        f.readline()
                        f.readline()
                        f.readline()
                        for line in f:
                            splitted_line = line.split("\t")
                            frequencies[splitted_line[1]][splitted_line[2]] = int(splitted_line[5])
                            if int(splitted_line[3]) - int(splitted_line[2]) > 2:
                                frequencies[splitted_line[1]][splitted_line[3]] = int(splitted_line[5])

    for root, dirs, files in os.walk(results_folder):
        for file in files:
            if file.endswith(".vcf.gz"):
                counter += 1
                sample_name = (file[:7])
                print(sample_name)
                if sample_name in QC_passed_samples:
                    for folder in folders:
                        to_output = []
                        header, cnvs = extract_cnvs(root + "/" + file)

                        header[-1] += "\t" + "genes_from_ERN" + "\t" + "Left/Right_breakpoint_frequency"
                        gene_list = gene_lists[folder]
                        for cnv in cnvs:
                            copy_of_cnv = [elem for elem in cnv]
                            if copy_of_cnv[6] == "MinQUAL":
                                continue
                            end_of_cnv = copy_of_cnv[1]
                            type_of_SV = copy_of_cnv[2].split(":")[0]
                            if copy_of_cnv[7].startswith("END="):
                                end_of_cnv = copy_of_cnv[7].split(";")[0][4:]
                            #if "[" in copy_of_cnv[7] or "]" in copy_of_cnv[7]:
                            #    end_of_cnv = -1
                            if not type_of_SV in number_of_SVs:
                                number_of_SVs[type_of_SV] = defaultdict(int)
                            if not type_of_SV in names_of_samples:
                                names_of_samples[type_of_SV] = defaultdict(set)
                            number_of_SVs[type_of_SV][(copy_of_cnv[0], copy_of_cnv[1], end_of_cnv)] += 1
                            names_of_samples[type_of_SV][(copy_of_cnv[0], copy_of_cnv[1], end_of_cnv)].add(sample_name)

                            cnv_coords = (copy_of_cnv[0], int(copy_of_cnv[1]), int(end_of_cnv))
                            freq_pair = "NA&NA"


                            intersection_list = []
                            for gene in gene_list[cnv_coords[0]]:
                                if are_intersected(cnv_coords, gene):
                                    intersection_list.append(gene[-1])
                            if len(intersection_list) > 0:
                                to_print = True
                                if annotate_frequencies:
                                    freq_pair = [0, 0]
                                    length_around = 20

                                    for i in range(int(copy_of_cnv[1]) - length_around, int(copy_of_cnv[1]) + length_around + 1):
                                        if str(i) in frequencies[copy_of_cnv[0]]:
                                            freq_pair[0] += frequencies[copy_of_cnv[0]][str(i)]
                                    if abs(int(end_of_cnv) - int(copy_of_cnv[1]) < 5):
                                        freq_pair[1] = freq_pair[0]
                                    else:
                                        for i in range(int(end_of_cnv) - length_around, int(end_of_cnv) + length_around + 1):
                                            if str(i) in frequencies[copy_of_cnv[0]]:
                                                freq_pair[1] += frequencies[copy_of_cnv[0]][str(i)]
                                    if freq_pair[0] > 20 or freq_pair[1] > 20:
                                        to_print = False
                                    freq_pair = "&".join(map(str, freq_pair))
                                copy_of_cnv.append(",".join(set(intersection_list)))
                                copy_of_cnv.append(freq_pair)
                                if to_print:
                                    to_output.append("\t".join(copy_of_cnv))
                        to_output = set(to_output)

                        if len(to_output) > 0 and len(to_output) < 10:
                            with open(path_to_output_folder + folder + "/" + sample_name + "_" + QC_passed_samples[sample_name] + "_exome_SVs.tsv", "w") as f:
                                f.write(header[-1] + "\n")
                                for elem in to_output:
                                    f.write(elem + "\n")
    if calculate_frequencies:
        for key in number_of_SVs:
            with open("/path/frequencies/" + key + ".seg", "w") as f:
                f.write("#type=GENE_EXPRESSION\n")
                f.write("#track graphtype=points name=\"" + key + "\" color=0,0,255 altColor=255,0,0 maxHeightPixels=80:80:80 viewLimits=0:0:" + str(counter) + " yLineMark=1 yLineOnOff=on\n")
                f.write("ID\tchr\tstart\tend\tsample\tvalue\n")
                for sv in names_of_samples[key]:
                    samples_with_sv = ";".join(names_of_samples[key][sv])
                    if sv[0].startswith("chrchr"):
                        print(samples_with_sv)
                    to_output_line = "\t".join((key, sv[0], sv[1], sv[2], samples_with_sv, str(len(names_of_samples[key][sv]))))
                    f.write(to_output_line + "\n")








def main():
    gene_list_ithaca = "/path/ERN-ITHACA.2021-06-23.3081genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv"
    gene_list_genturis = "/path/ERN-GENTURIS.2021-06-03.230genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv"
    gene_list_nmd = "/path/MuscleGeneTable.2021-05-26.611genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv"
    gene_list_rnd = "/path/ERN-RND.2021-07-13.1820genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv"

    folders = ["ithaca", "genturis", "nmd", "rmd"]
    path_to_output_folder = "/path/SVs/analysis/final_results/"
    QC_passed_samples = "/path/SVs/merged_1st_2nd.tsv"


    list_of_gene_lists = [gene_list_ithaca, gene_list_genturis, gene_list_nmd, gene_list_rnd]
    results_folder = "/path/MantaResults/"
    calculate_frequencies = True
    annotate_frequencies = False
    intersect_samples(list_of_gene_lists, results_folder, folders, path_to_output_folder, QC_passed_samples, calculate_frequencies, annotate_frequencies)

main()
