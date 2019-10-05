# -*- coding: utf-8 -*-

"""longdist.longdist: provides entry point main()."""

import os
import re
import csv
import numpy as npy
import configparser
import matplotlib.pyplot as plt

from argparse import ArgumentParser
from .sequence_attributes import SequenceAttributes
from sklearn import svm
from sklearn.model_selection import cross_val_score
from .pca_attributes import PCAAttributes
from multiprocessing import Pool
from math import floor
from sklearn import metrics
from sklearn.externals import joblib
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

__version__ = "1.0.3"


def main(args):
    '''
    parser = ArgumentParser(
        description='longdist: Method implementation for long ncRNAs and PCT distinction. This application can create '
                    'and use models base on the method by Schneider et al (2017).')
    parser.add_argument('--citation', action='store_true', help='Prints bibtex citation.')
    parser.add_argument('--version', action='store_true', help='Prints version number.')
    group = parser.add_argument_group("Method Paramenters")
    group.add_argument('--longs', metavar='<longs1.fa longs2.fa ...>', dest='longs', nargs='+',
                       help='List of fasta files containing long non-coding RNAs. The files should the separated by '
                            'spaces and should be ordered by species in same order as the pct file list. This '
                            'argument is required.')
    group.add_argument('--pcts', metavar='<pcts1.fa pcts2.fa ...>', dest='pcts', nargs='+',
                       help='List of fasta files containing protein coding transcripts. The files should the '
                            'separated by spaces and should be ordered by species in same order as the lcnRNA file '
                            'list. This argument is required.')
    group.add_argument('--input', nargs=1, metavar='<input.fa>', dest='input',
                       help='Fasta file containing transcripts to predict with the model')

    group.add_argument('--kmers', nargs=1, metavar='<50>', default=[50], type=int, dest='kmers',
                       help='Number of nucleotide pattern frequencies to consider in the model. Default is 50.')

    group.add_argument('--orf', nargs=1, metavar='<1>', default=[1], type=int, dest='orf',
                       help='The orf feature to be used by the model. Default is 1. Possible values are: 0 - No '
                            'orf feature; 1 - First ORF relative length; 2 - Longest ORF relative length')

    group.add_argument('--ratio', nargs=1, metavar='<0.75>', default=[0.75], type=float, dest='fraction',
                       help='The ratio of whole dataset that should be used for training. Default is 0.75.')

    group.add_argument('--size', nargs=1, metavar='<200>', default=[200], type=int, dest='size',
                       help='Minimum sequence size to consider. Default is 200.')
    group.add_argument('--cv', nargs=1, metavar='<10>', default=[10], type=int, dest='cross_validation',
                       help='Number of folds in cross-validation. Default is 10.')

    group.add_argument('--log2c', nargs=1, metavar='<-5,15,2>', default=["-5,15,2"],
                       help='Set the range of c to 2^{begin,...,begin+k*step,...,end}. Default is -5,15,2.')

    group.add_argument('--log2g', nargs=1, metavar='<3,-15,-2>', default=["3,-15,-2"],
                       help='Set the range of g to 2^{begin,...,begin+k*step,...,end}. Default is 3,-15,-2.')

    group.add_argument('--processes', nargs=1, metavar='<1>', default=[1], type=int,
                       help='Number of parallel processes for parameters search. Default is 1.')

    group.add_argument('--out_roc', nargs=1, metavar='<"lncRNA file"x"PCT file"x"kmers"_roc.eps>', dest='roc_file',
                       help='Name of the output file for the roc Curve. Default is roc.eps.')

    group.add_argument('--out_csv', nargs=1, metavar='<"lncRNA file"x"PCT file"x"kmers".csv>', dest='csv_file',
                       help='Name of the output CSV file containing the results. Default is a name built from the names '
                            'of both fasta files.')

    group.add_argument('--out_model', nargs=1, metavar='<"lncRNA file"x"PCT file"x"kmers".plk>',
                       dest='model_file',
                       help='Name of the output file containing the SVM Model. Default is a name built from the names '
                            'of both fasta files.')

    group.add_argument('--predict', action="store_true",
                       help='Just use a predefined model to distinguish long ncRNAs and PCTs in the input fasta file')

    group.add_argument('--model_config', nargs=1, metavar='<"lncRNA file"x"PCT file"x"kmers".plk>',
                       dest='model_config',
                       help='The file name containing the model configuration properties for prediction.')

    group.add_argument('--out', nargs=1, metavar='<"Input File".csv>',
                       dest='output',
                       help='The output CSV file for the result of the distinction made in the input file.')

    group.add_argument('--purge', action="store_true",
                       help='Purge all intermediate files. Intermediate files have ".longdist." in their names.'
                            ' All intermediate files are used to accelerate consecutive runs of the method. '
                            'Don\'t purge them if you want to run this method a second time with the same data')

    args = parser.parse_args()
    if args.citation:
        print("""@artile {Schneider:2017, title={A Support Vector Machine based method to distinguish long non-coding 
        RNAs from protein coding transcripts}, author={Schneider, Hugo and Raiol, Tainá and Brígido, Marcelo and 
        Walter, Maria E. M. T. and Stadler, Peter }, year={2017} }""")
    elif args.version:
        print(parser.description)
        print("Version: %s" % __version__)
    if args.predict:
        if args.input and args.species:
            predict(args)
        else:
            print("--input and --species are required for prediction")
    elif args.longs and args.pcts:
        print(parser.description)
        create_model(args)
    else:
        parser.print_usage()
	'''
    predict(args);


def predict(args):
    '''
    if not os.path.isfile(args.model_config[0]):
        print("Invalid model configuration file.")
        exit(1)
    '''
    config = configparser.ConfigParser()
    if args.species == 'Human':
	    args.model_config = 'GRCh38_GRCm38_firstOrf.plk.conf'
	    config.read(r'models/longdist/models/GRCh38_GRCm38_firstOrf.plk.conf')
	    clf = joblib.load(r'models/longdist/models/GRCh38_GRCm38_firstOrf.plk')
    elif args.species == 'Mouse':
	    args.model_config = 'GRCm38_GRCz10_firstOrf.plk.conf'
	    config.read(r'models/longdist/models/GRCm38_GRCz10_firstOrf.plk.conf')
	    clf = joblib.load(r'models/longdist/models/GRCm38_GRCz10_firstOrf.plk')
    kmers = config['MODEL']['attributes']
    kmers = re.sub("'\s+'", "', '", kmers) # wrong list format
    kmers = eval(kmers)
    fasta_input = SequenceAttributes(input_file=args.fasta, size=200, clazz=-1, use_intermediate_file=False)
    fasta_input.process(kmers)

    #clf = joblib.load(os.path.join(os.path.split(args.model_config)[0], config['MODEL']['model']))

    x = fasta_input.data[npy.array(kmers)]
    x = npy.array([list(l) for l in x])

    probabilities = clf.predict_proba(x)

    #csv_file = args.output[0] if args.output else "%s.csv" % args.input[0]
    csv_file = str(args.outfile)

    dump_result_csv(fasta_input.data["id"], probabilities[:, 1], probabilities[:, 0], csv_file)

    if args.purge:
        purge([fasta_input.intermediate_file()])


def features(args, kmers):
    if args.orf:
        if args.orf[0] == 0:
            return npy.array(kmers)
        elif args.orf[0] == 1:
            return npy.array(["fp"] + kmers)
        elif args.orf[0] == 2:
            return npy.array(["lp"] + kmers)
        else:
            print("Invalid value parameter ORF.")
            exit(1)
    else:
        return npy.array(["fp"] + kmers)


def create_model(args):
    longs = []
    pcts = []
    print("Processing fasta files. This could take some minutes... (if you don't have some intermediate files)")

    training = None
    testing = None
    pca = None

    if args.longs == args.pcts:
        print("Is not allowed to input the same list of files for lncRNAs and PCTs.")
        exit(1)

    for (long, pct) in zip(args.longs, args.pcts):
        l = SequenceAttributes(input_file=long, size=args.size[0], clazz=1)
        p = SequenceAttributes(input_file=pct, size=args.size[0], clazz=0)

        longs.append(l)
        pcts.append(p)

        print("Processing long non-coding RNA fasta file '%s'..." % long)
        l.process()
        print("Processing proteing coding transcripts fasta file '%s'..." % pct)
        p.process()

        min_size = min([len(l.data), len(p.data)])

        if pca is None:
            pca = npy.hstack((l.data, p.data))
        else:
            pca = npy.hstack((pca, npy.hstack((l.data, p.data))))

        longs_data_training, longs_data_testing = section(l.data, min_size, args.fraction[0])
        pcts_data_training, pcts_data_testing = section(p.data, min_size, args.fraction[0])

        if training is None:
            training = npy.hstack((longs_data_training, pcts_data_training))
        else:
            training = npy.hstack((training, npy.hstack((longs_data_training, pcts_data_training))))

        if testing is None:
            testing = npy.hstack((longs_data_testing, pcts_data_testing))
        else:
            testing = npy.hstack((testing, npy.hstack((longs_data_testing, pcts_data_testing))))

    pca = PCAAttributes(data=pca, patterns=SequenceAttributes.ALL_PATTERNS)
    kmers = pca.attributes(args.kmers[0])

    f = features(args, kmers)

    print("Selected features are: %s" % f)

    labels = training["class"]
    attributes = npy.array([list(l) for l in training[f]])

    testing_labels = testing["class"]
    testing_attributes = npy.array([list(l) for l in testing[f]])

    base_name = build_base_name(args.longs, args.pcts, args.kmers[0], args.orf[0])
    grid_file_name = "%s.longdist.npy" % base_name

    model_file = args.model_file[0] if args.model_file else "%s.plk" % base_name
    model_config_file = "%s.conf" % model_file
    csv_file = args.csv_file[0] if args.csv_file else "%s.csv" % base_name
    roc_file = args.roc_file[0] if args.roc_file else "%s_roc.eps" % base_name

    if os.path.exists(model_file):
        print("Using pre-built model ...")
        clf = joblib.load(model_file)
    else:
        c, gamma = svm_model_selection(attributes, labels, args.cross_validation[0], args.log2c, args.log2g,
                                       args.processes[0], grid_file_name)
        print("Building the model ...")
        clf = create_classifier(c=c, gamma=gamma)
        clf.fit(attributes, labels)
        joblib.dump(clf, model_file)

        config = configparser.ConfigParser()
        config['MODEL'] = {
            'desc': "Model built with lncRNA data from '%s' and PCT data from '%s'" % (
                os.path.basename(args.longs[0]), os.path.basename(args.pcts[0])),
            'attributes': npy.asarray(f),
            'model': os.path.relpath(model_file, os.path.split(model_config_file)[0])
        }

        with open(model_config_file, 'w') as config_file:
            config.write(config_file)
            args.model_config = [model_config_file]

    probabilities = clf.predict_proba(testing_attributes)
    long_probabilities = probabilities[:, 1]

    false_positive_rate, true_positive_rate, _ = metrics.roc_curve(testing_labels, long_probabilities)
    auc = metrics.auc(false_positive_rate, true_positive_rate)
    accuracy, sensitivity, specificity = accuracy_sensitivity_specificity(testing_labels, long_probabilities)

    dump_result_csv(testing["id"], long_probabilities, probabilities[:, 0], csv_file)
    roc(false_positive_rate, true_positive_rate, "AUC: %.2f%%" % (auc * 100),
        "%s x %s\nAccuracy: %.2f%% | Sensitivity: %.2f%% | Specificity: %.2f%%" % (
            os.path.basename(args.longs[0]), os.path.basename(args.pcts[0]), 100 * accuracy, 100 * sensitivity,
            100 * specificity), roc_file)

    if args.purge:
        purge([grid_file_name] + [file.intermediate_file() for file in longs] + [file.intermediate_file() for file in
                                                                                 pcts])

    if args.input:
        predict(args)


def purge(files):
    for f in files:
        os.remove(f)


def roc(false_positive_rate, true_positive_rate, label, title, file_name):
    fig, ax = plt.subplots()
    axins = zoomed_inset_axes(ax, 3.2, loc=7)

    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])

    ax.plot(false_positive_rate, true_positive_rate, label=label)
    axins.plot(false_positive_rate, true_positive_rate, label=label)

    axins.set_xlim(0.0, 0.1)  # apply the x-limits
    axins.set_ylim(0.9, 1)

    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")

    ax.plot([0, 1], [0, 1], 'k--')

    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title(title)
    ax.legend(loc=4)
    plt.savefig(file_name, format='eps')


def accuracy_sensitivity_specificity(labels, probabilities):
    pred = npy.copy(probabilities)
    pred[pred < 0.5] = 0
    pred[pred >= 0.5] = 1

    accuracy = metrics.accuracy_score(labels, pred)
    confusion_matrix = metrics.confusion_matrix(labels, pred)
    tp = confusion_matrix[1, 1]
    tn = confusion_matrix[0, 0]
    fp = confusion_matrix[0, 1]
    fn = confusion_matrix[1, 0]
    sensitivity = float(tp) / (float(fn + tp) if float(fn + tp) > 0 else -1)
    specificity = float(tn) / (float(tn + fp) if float(tn + fp) > 0 else -1)

    return accuracy, sensitivity, specificity


def dump_result_csv(ids, long_probabilities, pct_probabilities, file):
    with open(file, 'w') as csvfile:
        fieldnames = ['sequence', 'pct %', 'lncRNA %']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for (ID, lp, pp) in zip(ids, long_probabilities, pct_probabilities):
            writer.writerow({'sequence': ID, 'pct %': pp, 'lncRNA %': lp})


def build_base_name(long_files, pct_files, kmers, orf):
    long_file_base = []
    for long_file in long_files:
        long_file_base.append('.'.join(os.path.basename(long_file).split(sep='.')[:-1]))

    long_file_base = "_x_".join(long_file_base)

    pct_file_base = []
    for pct_file in pct_files:
        pct_file_base.append('.'.join(os.path.basename(pct_file).split(sep='.')[:-1]))

    pct_file_base = "_x_".join(pct_file_base)

    return os.path.join(os.path.split(long_files[0])[0],
                        "%s_x_%s_%d_orf%d" % (long_file_base, pct_file_base, kmers, orf))


def svm_model_selection(attributes, labels, folds, log2c, log2g, processes, file_name):
    print("Starting the SVM parameter search ...")
    c_begin, c_end, c_step = map(int, log2c[0].split(','))
    g_begin, g_end, g_step = map(int, log2g[0].split(','))

    pool = Pool(processes)

    if os.path.exists(file_name):
        results = npy.load(file_name).tolist()
    else:
        results = []

    def callback(result):
        results.append(result)

        print("C=%.13f, Gamma=%.13f: Accuracy=%.13f" % (
            result[0], result[1], result[2]))
        npy.save(file_name, npy.array(results))

    for log2c in range(c_begin, c_end, c_step):
        for log2g in range(g_begin, g_end, g_step):
            c, gamma = 2 ** log2c, 2 ** log2g
            if len(results) > 0:
                n_array = npy.array(results)
                index = npy.where(npy.all(n_array[:, :2] == npy.array([c, gamma]), axis=1))
                if len(index) > 0 and len(index[0]):
                    print("C=%.13f, Gamma=%.13f: Accuracy=%.13f (Restored from intermediate file)" % (
                        c, gamma, n_array[index[0], 2]))
                    continue

            pool.apply_async(cross_validation, args=(c, gamma, attributes, labels, folds),
                             callback=callback)

    pool.close()
    pool.join()

    results = npy.array(results)
    npy.save(file_name, results)

    best_models = results[results[:, 2] == npy.amax(results[:, 2])]

    best_models = best_models[best_models[:, 0] == npy.amin(best_models[:, 0])]
    [c, gamma, _] = best_models[0]

    return c, gamma


def create_classifier(c, gamma, verbose=False, shrinking=True, probability=True):
    return svm.SVC(kernel='rbf', C=c, gamma=gamma, decision_function_shape='ovr', max_iter=-1, tol=0.001,
                   verbose=verbose, shrinking=shrinking, probability=probability)


def cross_validation(c, gamma, attributes, labels, folds):
    clf = create_classifier(c=c, gamma=gamma, shrinking=False, probability=False)
    scores = cross_val_score(clf, attributes, labels, cv=folds)

    return [c, gamma, npy.max(scores)]


def section(data, size, fraction):
    if len(data) == size:
        remaining_data = data
    else:
        idx = npy.random.randint(len(data), size=size)
        remaining_data = data[idx]

    idx = npy.random.randint(size, size=int(floor(size * fraction)))
    mask = npy.ones(size, npy.bool)
    mask[idx] = 0

    return remaining_data[idx], remaining_data[mask]
