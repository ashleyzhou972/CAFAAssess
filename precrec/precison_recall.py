
def precision(tp,fp):
    return tp / (tp+fp)

def recall(tp,fn):
    return tp/(tp+fn)

def f1(prec, rec):
    return 2 * prec*rec/(prec+rec)

# Threshold dependent precision-recall


