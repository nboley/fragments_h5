from fragments_h5.fragment import Fragment


def test_fragment_eq():
    f1 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=None)
    f2 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=None)
    assert f1 == f2

    f1 = Fragment(chrom="chr1", start=629257, stop=629421)
    f2 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=None)
    assert f1 != f2

    f1 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=None)
    f2 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13)
    assert f1 == f2

    f1 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=None)
    f2 = Fragment(chrom="chr1", start=629257, stop=629421, mapq1=27, mapq2=13, gc=0.5)
    assert f1 != f2
