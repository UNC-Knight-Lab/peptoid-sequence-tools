import numpy as np
import pandas as pd
from itertools import groupby

################################################################
###################### Global Variables ########################
################################################################

G = 57.02146  # mass of glycine
B = 113.08406  # mass of N-butylglycine
Na = 22.989218  # mass of sodium
L = 855.55285  # linker mass for all other sequences.
d = 300  # discrimination value in ppm
ppm = [[]]  # allowed ppm error
N = 20  # number of monomers


################################################################
######################## Main function #########################
################################################################

def seq_match():
    print("Start evaluation.")
    peaks = data_extraction()
    y_ions = sequence_matching(peaks)
    export_to_excel(y_ions)
    print("Evaluation complete.")


################################################################
##################### Data Import #######################
################################################################

def data_extraction():
    print("Collecting data from folder...")
    peaks_list = []

    peaks_list.append(np.loadtxt("/Users/suprajachittari/Documents/python_scripts/Erin/J21_2.txt"))
    peaks_list.append(np.loadtxt("/Users/suprajachittari/Documents/python_scripts/Erin/J21_3.txt"))
    peaks_list.append(np.loadtxt("/Users/suprajachittari/Documents/python_scripts/Erin/J21_4.txt"))

    peaks = np.concatenate(peaks_list)
    return peaks


################################################################
##################### Matching Functions #######################
################################################################

def match_two(curr_sum, peaks, i):
    two_diffs = []
    species = [G + G + Na, G + B + Na, B + B + Na]
    species_label1 = ['G+Na', 'G or B', 'B+Na']
    species_label2 = ['G+Na', 'G or B', 'B+Na']

    for f in range(len(species)):
        for p in peaks:
            test_diff = p - curr_sum - species[f]
            if abs(test_diff) <= d * ((curr_sum + species[f]) / 1E6):
                match1 = {
                    "position": i + 1,
                    "species": species_label1[f],
                    "difference": test_diff,
                    "target peak": 0,
                    "picked peak": 0,
                    "no adduct": 0
                }

                match2 = {
                    "position": i + 2,
                    "species": species_label2[f],
                    "difference": test_diff,
                    "target peak": curr_sum + species[f],
                    "picked peak": p,
                    "no adduct": curr_sum + species[f] - Na
                }

                two_diffs.append(match1)
                two_diffs.append(match2)
    return two_diffs


def match_one(curr_sum, peaks, i):
    diffs = []

    species_label = ['G+Na', 'B+Na']
    species = [G + Na, B + Na]

    for f in range(len(species)):
        for p in peaks:
            test_diff = p - curr_sum - species[f]
            if abs(test_diff) <= d * ((curr_sum + species[f]) / 1E6):
                match = {
                    "position": i + 1,
                    "species": species_label[f],
                    "difference": test_diff,
                    "target peak": curr_sum + species[f],
                    "picked peak": p,
                    "no adduct": curr_sum + species[f] - Na
                }
                diffs.append(match)
    return diffs


################################################################
##################### Sequence matching function ###############
################################################################

def sequence_matching(peaks):
    print("Beginning sequence matching.")
    Y_ions = [[]]  # y-ions list

    print("Appending linker mass to start of chain.")
    # Append linker mass
    Y_ions[0].append({
        "position": 0,
        "species": "linker",
        "difference": 0,
        "target peak": L + Na,
        "picked peak": 0,
        "no adduct": L
    })

    for i in range(N):  # iterate over all positions in chain
        tot = len(Y_ions)
        new_lists = []
        remove_nums = []
        for it in range(tot):
            if Y_ions[it][-1]["position"] == i:
                print("Evaluating position", i + 1, "for sequence ", it + 1)

                # print("Evaluating sequence ", it+1)
                curr_sum = Y_ions[it][-1]['no adduct']

                # compare each predicted Y ion + adducts to individual peak - guessing single matching
                print("Attempting single residue matching.")
                diffs = match_one(curr_sum, peaks, i)

                # filter diffs list
                diffs_sorted = sorted(diffs, key=lambda d: d['difference'])
                diffs_unique = [next(d) for _, d in groupby(diffs_sorted, key=lambda _d: _d['species'])]

                # guess two
                if len(diffs_unique) == 0:
                    diffs_mult = match_two(curr_sum, peaks, i)

                    if len(diffs_mult) != 0:
                        print("Single residuce matching failed. Pair residuce matching was success.")
                        Y_ions[it].append(diffs_mult[0])
                        Y_ions[it].append(diffs_mult[1])
                    else:
                        print("Single residue matching failed. Pair residue matching failed. The branch is dead.")

                # offer user input and branch to multiple sequence options
                elif len(diffs_unique) > 1:
                    G_candidates = []
                    B_candidates = []

                    print("More than one sequence match was found. The following are the options:")

                    for jj in range(len(diffs)):
                        print("Species", diffs[jj]["species"], "has difference", round(diffs[jj]["difference"], 4),
                              "picked peak:", round(diffs[jj]["picked peak"], 4), "target peak:",
                              round(diffs[jj]["target peak"], 4))
                        if diffs[jj]["species"] == "G+Na":
                            G_candidates.append(diffs[jj])
                        elif diffs[jj]["species"] == "B+Na":
                            B_candidates.append(diffs[jj])

                    while True:
                        try:
                            picked = input(
                                "Enter 'G' to choose glycine, 'B' to choose N-butyl, or 'G B' to choose to branch:")
                            if all(x.isalpha() or x.isspace() for x in picked):
                                picked = list(map(str, picked.split(' ')))
                                break
                            else:
                                raise TypeError
                        except TypeError:
                            print("Sorry. Please input either 'G', 'B', or 'G B'.")
                            continue
                        except EOFError:
                            print("Sorry. Please input either 'G', 'B', or 'G B'.")
                            continue

                    minG = min(G_candidates, key=lambda x: x['difference'])
                    minB = min(B_candidates, key=lambda x: x['difference'])

                    copy_y = Y_ions[it].copy()

                    for kk in range(len(picked)):
                        if picked[kk] == "G":  # assign G as the next position with the smallest ppm difference
                            new_lists.append(copy_y + [minG])
                            print("You picked ", picked[kk])
                        else:  # assign B as the next position with the smallest ppm difference
                            new_lists.append(copy_y + [minB])
                            print("You picked ", picked[kk])

                    remove_nums.append(it)

                # assign the only choice
                else:
                    print("Single matching success.")
                    Y_ions[it].append(diffs_unique[0])

        # remove duplicate sequences
        for j in sorted(remove_nums, reverse=True):
            Y_ions.remove(Y_ions[j])

        Y_ions += new_lists

        print("Current number of sequences:", len(Y_ions))

    print("Sequence matching complete.")
    return Y_ions


################################################################
######################## Excel export ##########################
################################################################

def export_to_excel(Y_ions):
    print("Exporting data to excel in directory.")
    writer = pd.ExcelWriter('seq_1.xlsx', engine='xlsxwriter')

    for i in range(len(Y_ions)):
        df = pd.DataFrame(Y_ions[i])
        df.to_excel(writer, sheet_name='Sheet{}'.format(i + 1))
    writer.save()