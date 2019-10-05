from Bio import SeqIO
import numpy as np
import re
import os


class SequenceAttributes:
    LONGDIST_NPY = "%s.longdist.npy"

    ALL_PATTERNS = [
        "fl",
        "fp",
        "ll",
        "lp",
        "aa",
        "ac",
        "ag",
        "at",
        "ca",
        "cc",
        "cg",
        "ct",
        "ga",
        "gc",
        "gg",
        "gt",
        "ta",
        "tc",
        "tg",
        "tt",
        "aaa",
        "aac",
        "aag",
        "aat",
        "aca",
        "acc",
        "acg",
        "act",
        "aga",
        "agc",
        "agg",
        "agt",
        "ata",
        "atc",
        "atg",
        "att",
        "caa",
        "cac",
        "cag",
        "cat",
        "cca",
        "ccc",
        "ccg",
        "cct",
        "cga",
        "cgc",
        "cgg",
        "cgt",
        "cta",
        "ctc",
        "ctg",
        "ctt",
        "gaa",
        "gac",
        "gag",
        "gat",
        "gca",
        "gcc",
        "gcg",
        "gct",
        "gga",
        "ggc",
        "ggg",
        "ggt",
        "gta",
        "gtc",
        "gtg",
        "gtt",
        "taa",
        "tac",
        "tag",
        "tat",
        "tca",
        "tcc",
        "tcg",
        "tct",
        "tga",
        "tgc",
        "tgg",
        "tgt",
        "tta",
        "ttc",
        "ttg",
        "ttt",
        "aaaa",
        "aaac",
        "aaag",
        "aaat",
        "aaca",
        "aacc",
        "aacg",
        "aact",
        "aaga",
        "aagc",
        "aagg",
        "aagt",
        "aata",
        "aatc",
        "aatg",
        "aatt",
        "acaa",
        "acac",
        "acag",
        "acat",
        "acca",
        "accc",
        "accg",
        "acct",
        "acga",
        "acgc",
        "acgg",
        "acgt",
        "acta",
        "actc",
        "actg",
        "actt",
        "agaa",
        "agac",
        "agag",
        "agat",
        "agca",
        "agcc",
        "agcg",
        "agct",
        "agga",
        "aggc",
        "aggg",
        "aggt",
        "agta",
        "agtc",
        "agtg",
        "agtt",
        "ataa",
        "atac",
        "atag",
        "atat",
        "atca",
        "atcc",
        "atcg",
        "atct",
        "atga",
        "atgc",
        "atgg",
        "atgt",
        "atta",
        "attc",
        "attg",
        "attt",
        "caaa",
        "caac",
        "caag",
        "caat",
        "caca",
        "cacc",
        "cacg",
        "cact",
        "caga",
        "cagc",
        "cagg",
        "cagt",
        "cata",
        "catc",
        "catg",
        "catt",
        "ccaa",
        "ccac",
        "ccag",
        "ccat",
        "ccca",
        "cccc",
        "cccg",
        "ccct",
        "ccga",
        "ccgc",
        "ccgg",
        "ccgt",
        "ccta",
        "cctc",
        "cctg",
        "cctt",
        "cgaa",
        "cgac",
        "cgag",
        "cgat",
        "cgca",
        "cgcc",
        "cgcg",
        "cgct",
        "cgga",
        "cggc",
        "cggg",
        "cggt",
        "cgta",
        "cgtc",
        "cgtg",
        "cgtt",
        "ctaa",
        "ctac",
        "ctag",
        "ctat",
        "ctca",
        "ctcc",
        "ctcg",
        "ctct",
        "ctga",
        "ctgc",
        "ctgg",
        "ctgt",
        "ctta",
        "cttc",
        "cttg",
        "cttt",
        "gaaa",
        "gaac",
        "gaag",
        "gaat",
        "gaca",
        "gacc",
        "gacg",
        "gact",
        "gaga",
        "gagc",
        "gagg",
        "gagt",
        "gata",
        "gatc",
        "gatg",
        "gatt",
        "gcaa",
        "gcac",
        "gcag",
        "gcat",
        "gcca",
        "gccc",
        "gccg",
        "gcct",
        "gcga",
        "gcgc",
        "gcgg",
        "gcgt",
        "gcta",
        "gctc",
        "gctg",
        "gctt",
        "ggaa",
        "ggac",
        "ggag",
        "ggat",
        "ggca",
        "ggcc",
        "ggcg",
        "ggct",
        "ggga",
        "gggc",
        "gggg",
        "gggt",
        "ggta",
        "ggtc",
        "ggtg",
        "ggtt",
        "gtaa",
        "gtac",
        "gtag",
        "gtat",
        "gtca",
        "gtcc",
        "gtcg",
        "gtct",
        "gtga",
        "gtgc",
        "gtgg",
        "gtgt",
        "gtta",
        "gttc",
        "gttg",
        "gttt",
        "taaa",
        "taac",
        "taag",
        "taat",
        "taca",
        "tacc",
        "tacg",
        "tact",
        "taga",
        "tagc",
        "tagg",
        "tagt",
        "tata",
        "tatc",
        "tatg",
        "tatt",
        "tcaa",
        "tcac",
        "tcag",
        "tcat",
        "tcca",
        "tccc",
        "tccg",
        "tcct",
        "tcga",
        "tcgc",
        "tcgg",
        "tcgt",
        "tcta",
        "tctc",
        "tctg",
        "tctt",
        "tgaa",
        "tgac",
        "tgag",
        "tgat",
        "tgca",
        "tgcc",
        "tgcg",
        "tgct",
        "tgga",
        "tggc",
        "tggg",
        "tggt",
        "tgta",
        "tgtc",
        "tgtg",
        "tgtt",
        "ttaa",
        "ttac",
        "ttag",
        "ttat",
        "ttca",
        "ttcc",
        "ttcg",
        "ttct",
        "ttga",
        "ttgc",
        "ttgg",
        "ttgt",
        "ttta",
        "tttc",
        "tttg",
        "tttt"
    ]

    DI_TRI_PATTERNS = [
        "fl",
        "fp",
        "ll",
        "lp",
        "aa",
        "ac",
        "ag",
        "at",
        "ca",
        "cc",
        "cg",
        "ct",
        "ga",
        "gc",
        "gg",
        "gt",
        "ta",
        "tc",
        "tg",
        "tt",
        "aaa",
        "aac",
        "aag",
        "aat",
        "aca",
        "acc",
        "acg",
        "act",
        "aga",
        "agc",
        "agg",
        "agt",
        "ata",
        "atc",
        "atg",
        "att",
        "caa",
        "cac",
        "cag",
        "cat",
        "cca",
        "ccc",
        "ccg",
        "cct",
        "cga",
        "cgc",
        "cgg",
        "cgt",
        "cta",
        "ctc",
        "ctg",
        "ctt",
        "gaa",
        "gac",
        "gag",
        "gat",
        "gca",
        "gcc",
        "gcg",
        "gct",
        "gga",
        "ggc",
        "ggg",
        "ggt",
        "gta",
        "gtc",
        "gtg",
        "gtt",
        "taa",
        "tac",
        "tag",
        "tat",
        "tca",
        "tcc",
        "tcg",
        "tct",
        "tga",
        "tgc",
        "tgg",
        "tgt",
        "tta",
        "ttc",
        "ttg",
        "ttt"
    ]

    def __init__(self, input_file, size, clazz, use_intermediate_file=True):
        self.fasta_file = input_file
        self.size = size
        self.clazz = clazz
        self.use_intermediate_file = use_intermediate_file

    def intermediate_file(self):
        return self.LONGDIST_NPY % self.fasta_file

    def process(self, patterns=ALL_PATTERNS):

        if self.use_intermediate_file and os.path.exists(self.LONGDIST_NPY % self.fasta_file):
            self.data = np.load(self.intermediate_file())
            return self.data

        data = []
        with open(self.fasta_file, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if len(record.seq) >= self.size:
                    data.append(self.attributes(record, self.clazz, patterns))

        dt = np.dtype(
            [("id", np.str_, 32), ("class", np.int_), ("length", np.int_)] + [(pattern, np.float64) for pattern in
                                                                              patterns])
        self.data = np.asarray(np.array(data, dtype=dt))
        if self.use_intermediate_file:
            self.dump()

        return self.data

    def attributes(self, record, type, patterns):
        first_orf_size = None
        longest_orf_size = None
        seq = str(record.seq)
        attributes = [record.id, type, len(seq)]
        for pattern in patterns:
            if pattern == "ll":
                if longest_orf_size == None:
                    longest_orf_size = self.longest_orf(seq)
                attributes.append(longest_orf_size)
            elif pattern == "lp":
                if longest_orf_size == None:
                    longest_orf_size = self.longest_orf(seq)
                attributes.append(float(longest_orf_size) / float(len(seq)))
            elif pattern == "fl":
                if first_orf_size == None:
                    first_orf_size = self.first_orf(seq)
                attributes.append(first_orf_size)
            elif pattern == "fp":
                if first_orf_size == None:
                    first_orf_size = self.first_orf(seq)
                attributes.append(float(first_orf_size) / float(len(seq)))
            else:
                attributes.append(self.count_pattern(pattern, seq))

        return tuple(attributes)

    def count_pattern(self, pattern, seq):
        length = len(seq)
        count = len([m.start() for m in re.finditer("(?=%s)" % pattern, seq, re.IGNORECASE)])
        total = 0
        attr_length = len(pattern)
        for j in range(0, len(pattern)):
            total += int((length - j) / attr_length)

        return float(count) / float(total)

    def first_orf(self, seq):
        index = re.search("atg", seq, re.IGNORECASE)
        if index == None:
            return 0
        else:
            return self.orf_size(index.start(), seq)

    def longest_orf(self, seq):
        sizes = []
        for index in [m.start() for m in re.finditer("atg", seq, re.IGNORECASE)]:
            sizes.append(self.orf_size(index, seq))
        if (len(sizes) == 0):
            return 0
        else:
            return max(sizes)

    def orf_size(self, start, seq):
        for end in [m.start() for m in re.finditer("taa|tga|tag", seq, re.IGNORECASE)]:
            if end < start:
                continue

            length = end - start + 3
            if length % 3 == 0:
                return length

        length = len(seq[start:])
        length -= len(seq[start:]) % 3
        return length

    def dump(self):
        np.save(self.intermediate_file(), self.data)
