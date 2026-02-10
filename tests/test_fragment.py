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


def test_fragment_length_and_midpoint():
    f = Fragment(chrom="chr1", start=100, stop=200)
    assert f.length == 100
    assert f.tlen == 100
    assert f.midpoint == 150

    # Test odd-length fragment (rounds down)
    f2 = Fragment(chrom="chr1", start=100, stop=201)
    assert f2.length == 101
    assert f2.midpoint == 150  # 100 + 101//2 = 150


def test_fragment_mapq_min():
    # Both mapqs set
    f = Fragment(chrom="chr1", start=0, stop=100, mapq1=30, mapq2=20)
    assert f.mapq12_min == 20

    # Only mapq1 set
    f = Fragment(chrom="chr1", start=0, stop=100, mapq1=30, mapq2=None)
    assert f.mapq12_min == 30

    # Only mapq2 set
    f = Fragment(chrom="chr1", start=0, stop=100, mapq1=None, mapq2=20)
    assert f.mapq12_min == 20

    # Neither set
    f = Fragment(chrom="chr1", start=0, stop=100)
    assert f.mapq12_min is None


def test_fragment_mapq_gte():
    f = Fragment(chrom="chr1", start=0, stop=100, mapq1=30, mapq2=20)
    assert f.mapq_gte(20) is True
    assert f.mapq_gte(21) is False

    # None mapqs always pass
    f = Fragment(chrom="chr1", start=0, stop=100)
    assert f.mapq_gte(60) is True


def test_fragment_replace():
    f1 = Fragment(chrom="chr1", start=100, stop=200, mapq1=30)
    f2 = f1.replace(start=150)
    assert f2.chrom == "chr1"
    assert f2.start == 150
    assert f2.stop == 200
    assert f2.mapq1 == 30


def test_fragment_strand_normalization():
    # Binary strand values should be normalized
    f1 = Fragment(chrom="chr1", start=0, stop=100, strand="+")
    f2 = Fragment(chrom="chr1", start=0, stop=100, strand=b"+")
    assert f1.strand == "+"
    assert f2.strand == "+"


def test_fragment_end_clipped_field():
    # fragment_end_clipped is optional; None when unknown
    f = Fragment(chrom="chr1", start=0, stop=100)
    assert f.fragment_end_clipped is None

    f = Fragment(chrom="chr1", start=0, stop=100, fragment_end_clipped=False)
    assert f.fragment_end_clipped is False

    f = Fragment(chrom="chr1", start=0, stop=100, fragment_end_clipped=True)
    assert f.fragment_end_clipped is True
    assert "fragment_end_clipped=True" in repr(f)
