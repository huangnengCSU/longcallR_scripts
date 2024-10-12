def find_motif_donor(seq, gt_truth):
    """
    :param seq: reference sequence
    :param gt_truth: truth gt sites from annotation
    :return:
    """
    true_samples = []
    negative_samples = []
    for i in range(len(seq)):
        if i - 2 < 0:
            continue
        if seq[i - 2:i] == "GT" or seq[i - 2:i] == "gt":
            ## candidate GT site
            if i + 1 in gt_truth:
                ## gt_truth is 1-based, i is 0-based
                try:
                    true_samples.append(seq[i - 2 - 100:i + 100])
                except IndexError:
                    continue
            else:
                try:
                    negative_samples.append(seq[i - 2 - 100:i + 100])
                except IndexError:
                    continue

        if i + 3 > len(seq):
            continue
        if seq[i + 1:i + 3] == "AC" or seq[i + 1:i + 3] == "ac":
            if i + 1 in gt_truth:
                try:
                    true_samples.append(seq[i - 99:i + 2 + 101])
                except IndexError:
                    continue
            else:
                try:
                    negative_samples.append(seq[i - 99:i + 2 + 101])
                except IndexError:
                    continue
    return true_samples, negative_samples
