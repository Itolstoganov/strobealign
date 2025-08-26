#!/usr/bin/env python3

import argparse
import pygit2
import os
import shutil
import subprocess

from collections import namedtuple


MC_REPO = "https://github.com/Itolstoganov/strobealign"

Run = namedtuple('Run', ['params', 'version', 'search_level', 'type', 'k'])


def compile(build_dir: str, versions: dict, nthreads: int):
    for version, branch in versions.items():
        bin_name = os.path.abspath(os.path.join(build_dir, "strobealign-{}".format(branch)))
        if os.path.exists(bin_name):
            continue
        repo_path = os.path.join(build_dir, version)
        version_repo = pygit2.clone_repository(MC_REPO, repo_path, checkout_branch=branch)
        os.chdir(repo_path)
        print(os.getcwd())
        cmake_cmd = " ".join(["cmake", "-B", "build", "-DENABLE_AVX=ON", "-DISAL=download"])
        subprocess.run(cmake_cmd, shell=True, check=True)
        print(os.getcwd(), "running make")
        make_cmd = " ".join(["make", "-j", str(nthreads), "-C", "build"])
        print(make_cmd)
        subprocess.run(make_cmd, shell=True)
        print(bin_name)
        shutil.copy2("build/strobealign", bin_name)
        assert(os.path.exists(bin_name))
        os.chdir("../../../")


def make_kmer_params(k: int, m: int, t: int):
    return "-ci1 -k{} -r -m{} -t{} -fm -cs10000000".format(k, m, t)


def make_parameters(versions, main_k: int, mem: int, threads: int):
    parameters = {}
    for k in [main_k, main_k // 3, main_k // 2, main_k * 2 // 3]:
        parameters["kmers\t{}".format(k)] = Run(params=make_kmer_params(k, mem, threads), version="kmc", search_level=1, type="kmer", k=k)
    w_min = (main_k // 2) // 2 * 3
    w_max = (main_k // 2) * 3
    parameters["randstrobes\t(2,{},{},{})".format(main_k // 2, w_min, w_max)] = Run(params="-k {} -s {} -l {} -u {}".format(main_k // 2, main_k // 2, w_min, w_max), \
        version="2_randstrobes", search_level=1, type="randstrobes", k=main_k // 2)
    parameters["multi-context\t(2,{},{},{})".format(main_k // 2, w_min, w_max)] = Run(params="-k {} -s {} -l {} -u {} --mcs".format(main_k // 2, main_k // 2, w_min, w_max), \
        version="2-mcs", search_level=2, type="mcs", k=main_k // 2)
    parameters["randstrobes\t(3,{},{},{})".format(main_k // 3, w_min, w_max)] = Run(params="-k {} -s {} -l {} -u {}".format(main_k // 3, main_k // 3, w_min, w_max), \
        version="3_randstrobes", search_level=1, type="randstrobes", k=main_k // 3)
    parameters["multi-context\t(3,{},{},{})".format(main_k // 3, w_min, w_max)] = Run(params="-k {} -s {} -l {} -u {} --mcs".format(main_k // 3, main_k // 3, w_min, w_max), \
        version="3-strobes-experimental", search_level=3, type="mcs", k=main_k // 3)
    return parameters


def parse_index_stats(index_outpath: str, search_level: int):
    results = {}
    with open(index_outpath, "r") as stats_handle:
        for line in stats_handle:
            if line.startswith("E-hits"):
                line_lst = line.strip().split(",")
                assert(len(line_lst) == search_level + 1)
                for _ in range(3 - search_level):
                    line_lst.append("-")
                if search_level == 3:
                    line_lst[2], line_lst[3] = line_lst[3], line_lst[2]
                print(line_lst)
                results["E-hits"] = ["{:.3f}".format(float(x)) if x != '-' else '-' for x in line_lst[1:]]
            if line.startswith("Unique"):
                line_lst = line.strip().split(",")
                assert(len(line_lst) == search_level + 1)
                for _ in range(3 - search_level):
                    line_lst.append("-")
                if search_level == 3:
                    line_lst[2], line_lst[3] = line_lst[3], line_lst[2]
                results["Unique"] = ["{:.3f}".format(float(x)) if x != '-' else '-' for x in line_lst[1:]]
            if line.startswith("Distinct"):
                line_lst = line.strip().split(",")
                assert(len(line_lst) == search_level + 1)
                for _ in range(3 - search_level):
                    line_lst.append("-")
                if search_level == 3:
                    line_lst[2], line_lst[3] = line_lst[3], line_lst[2]
                results["Distinct"] = ["{:.1f}".format(float(x) / 1000000.0) if x != '-' else '-' for x in line_lst[1:]]
    return results


def parse_kmc_stats(index_path):
    num_unique = 0
    num_distinct = 0
    sum_occ_sq = 0
    sum_occ = 0
    with open(index_path, "r") as in_handle:
        for line in in_handle:
            occ = int(line.strip().split()[0])
            num_seeds = int(line.strip().split()[1])
            sum_occ += occ * num_seeds
            sum_occ_sq += occ * occ * num_seeds
            if occ == 1:
                num_unique = num_seeds
            num_distinct += num_seeds
    results = {}
    results["E-hits"] = ["{:.3f}".format(float(sum_occ_sq) / float(sum_occ)), '-', '-']
    results["Unique"] = ["{:.3f}".format(float(num_unique) / float(num_distinct)), '-', '-']
    results["Distinct"] = ["{:.1f}".format(float(num_distinct) / 1000000.0), '-', '-']
    return results


def run_index(output_dir: str, build_dir: str, versions: dict, k_values: list[int], reference_path: str, left_path: str, right_path: str, num_threads: int, mem: int):
    table_results = {}
    index = 0
    read_length = 500
    for k in k_values:
        parameters = make_parameters(versions, k, mem, num_threads)
        for description, run in parameters.items():
            if run.type != "kmer":
                bin_path = os.path.join(build_dir, "strobealign-{}".format(run.version))
                print(bin_path)
                run_id = "run_" + run.params.replace(" ", "_").replace("-", "") + "_" + run.version.replace("-", "_")
                index_outpath = os.path.join(output_dir, run_id + ".stats")
                if not os.path.exists(index_outpath):
                    index_cmd = " ".join([bin_path, reference_path, left_path, right_path, "--index-statistics={}".format(index_outpath), "-t {}".format(num_threads), \
                                          "-i", "-r {}".format(read_length), run.params])
                    print(index_cmd)
                    subprocess.run(index_cmd, shell=True)
            if run.type == "kmer":
                run_id = "run_kmc_{}".format(run.k)
                max_occ = 100000
                index_outpath = os.path.join(output_dir, run_id + ".kmc.txt")
                tmpdir = os.path.join(output_dir, run_id + "_tmp")
                if not os.path.exists(index_outpath):
                    if not os.path.exists(tmpdir):
                        os.mkdir(tmpdir)
                    kmc_cmd = " ".join(["kmc", run.params, reference_path, os.path.join(output_dir, run_id), os.path.join(output_dir, run_id + "_tmp")])
                    dump_cmd = "kmc_tools transform {} histogram {} -cx{}".format(os.path.join(output_dir, run_id), index_outpath, max_occ)
                    subprocess.run(kmc_cmd, shell=True)
                    subprocess.run(dump_cmd, shell=True)

            run_results = parse_index_stats(index_outpath, run.search_level) if run.type != "kmer" else parse_kmc_stats(index_outpath)
            print(description, run_results)
            table_results[description] = run_results
            for i in range(run.search_level):
                # print(run_results)
                # print(run_results["Unique"][i], run_results["E-hits"][i])
                search_description = ""
                run_k = k
                if run.search_level > 1:
                    if i == 0:
                        search_description = " (full seed)"
                    elif i == 1:
                        search_description = " (first seed)"
                        run_k = k // run.search_level
                    else:
                        search_description = " (first and second seed)"
                        run_k = 2 * k // run.search_level
                index += 1
    return table_results


def print_results_table(uniqueness_results: dict, outpath: str):
    dataset_name = "CHM13"
    table_string = "\\begin{table}[]\n"
    metrics = ["\% Unique", "E-hits", "\\#Distinct (millions)"]
    num_columns = len(metrics) * 3 + 2
    column_string = "{" + "".join(["l" for _ in range(num_columns)]) + "}"
    table_string += "\\begin{tabular}" + column_string + "\n"
    table_string += " &  & " + "\\multicolumn{{{}}}{{c}}{{{}}}".format(num_columns - 2, dataset_name) + "\\\\ \\cline{{{}-{}}}\n".format(3, num_columns)
    table_string += " &  "
    for metric in metrics:
        table_string += "& \\multicolumn{{{}}}{{c}}{{{}}} ".format(3, metric)
    table_string += "\\\\ \\cline{{{}-{}}}\n".format(3, num_columns)
    table_string += " & "
    level_list = ["Full", "1st", "1st \& 2nd"]
    for _ in range(len(metrics)):
        table_string += " & "
        table_string += " & ".join(level_list)
    table_string += " \\\\ \\hline\n"
    # Generate result
    for protocol, results in uniqueness_results.items():
        table_string += protocol.replace("\t"," & ")
        table_string += " & "
        table_string += " & ".join(results["Unique"])
        table_string += " & "
        table_string += " & ".join(results["E-hits"])
        table_string += " & "
        table_string += " & ".join(results["Distinct"])
        table_string += " \\\\\n"
    table_string += "\\end{tabular}\n"
    table_string += "\\caption{}\n"
    table_string += "\\end{table}\n"
    with open(outpath, "w") as outhandle:
        outhandle.write(table_string)


def createparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', "-r", help="Reference genome")
    parser.add_argument('--threads', "-t", help="Number of threads", type=int, default=4)
    parser.add_argument('--mem', '-m', help="Max RAM (GB)", type=int, default=8)
    parser.add_argument('--k', '-k', help="main k value", type=int, default=60)
    parser.add_argument('--output', '-o', help="Output directory")
    return parser


parser = createparser()
args = parser.parse_args()

output_dir = args.output

if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.mkdir(output_dir)

build_dir = os.path.join(output_dir, "bin")
nthreads = args.threads

versions = {"3_randstrobes": "3_randstrobes", "2_mcs": "2-mcs", "2_randstrobes": "2_randstrobes", "3_mcs": "3-strobes-experimental"}

compile(build_dir, versions, nthreads)

k_values = [args.k]
left_path = os.path.join(args.output, "empty1.fq")
right_path = os.path.join(args.output, "empty2.fq")
open(left_path, 'a').close()
open(right_path, 'a').close()
table_results = run_index(output_dir, build_dir, versions, k_values, args.ref, left_path, right_path, nthreads, args.mem)

outpath = os.path.join(output_dir, "stats.tex")
print_results_table(table_results, outpath)