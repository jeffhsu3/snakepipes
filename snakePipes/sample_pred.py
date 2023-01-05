""" Functions for sample prediction 
"""
import re
import glob
import os

import fuzzywuzzy as fuzz

from snakePipes.common_functions import get_sample_names_bam, write_configfile


def predict_chip_dict(wdir, input_pattern_str, bamExt, fromBAM=None):
    """
    Predict a chip_dict from set of bam files
    ChIP input/control samples are identified from input_pattern (default: 'input')
    for each sample then the best input sample (by fuzzywuzzy score) is selected
    chip_dict is written as yaml to workflow workingdir
    predicts whether a sample is broad or narrow based on histone mark pattern
    """
    pat = "|".join(re.split(',| |\\||;', input_pattern_str))
    input_pat = r".*(" + pat + ")"
    clean_pat = r"" + pat + ""
    pat1 = re.compile(clean_pat, re.IGNORECASE)

    if fromBAM:
        infiles = sorted(glob.glob(os.path.join(fromBAM, '*' + bamExt)))
    else:
        infiles = sorted(glob.glob(os.path.join(wdir, 'filtered_bam/', '*.bam')))
    samples = get_sample_names_bam(infiles, bamExt)

    chip_dict_pred = {}
    chip_dict_pred["chip_dict"] = {}
    print("---------------------------------------------------------------------------------------")
    print("Predict Chip-seq sample configuration")
    print("---------------------------------------------------------------------------------------")
    print("\nSearch for Input/control samples...")

    input_samples = set([])
    for i in samples:
        if re.match(input_pat, i, re.IGNORECASE):
            print("...found: ", i)
            input_samples.add(i)

    print("\nTry to find corresponding ChIP samples...")

    for i in samples:
        if i in input_samples:
            continue

        print("\n sample: ", i,)
        matches_sim = {}
        for j in input_samples:
            c_clean = pat1.sub("", j)
            sim1 = fuzz.ratio(c_clean, i) + fuzz.partial_ratio(c_clean, i) + fuzz.token_sort_ratio(c_clean, i) + fuzz.token_set_ratio(c_clean, i)
            matches_sim[j] = sim1 / 4

        sim = 0
        final_matches = set([])
        for key, value in sorted(matches_sim.items(), key=lambda k: (k[1], k[0]), reverse=True):
            if value >= sim:
                final_matches.add(key)
                print("   top matching input sample by score: %s = %s" % (key, value))
                sim = value

        tmp = ':'.join(list(final_matches))

        if len(final_matches) > 1:
            tmp = "__PLEASE_SELECT_ONLY_ONE_CONTROL__:" + tmp
        elif len(final_matches) == 0:
            print("No control sample found!")

        chip_dict_pred["chip_dict"][i] = {}
        chip_dict_pred["chip_dict"][i]['control'] = tmp
        if re.match(".*(H3K4me1|H3K36me3|H3K9me3|H3K27me3).*", i, re.IGNORECASE):
            chip_dict_pred["chip_dict"][i]['broad'] = True
        else:
            chip_dict_pred["chip_dict"][i]['broad'] = False

    outfile = os.path.join(wdir, "chip_seq_sample_config.PREDICTED.yaml")
    write_configfile(outfile, chip_dict_pred)
    print("---------------------------------------------------------------------------------------")
    print("Chip-seq sample configuration is written to file ", outfile)
    print("Please check and modify this file - this is just a guess! Then run the workflow with it.")
    print("---------------------------------------------------------------------------------------")

