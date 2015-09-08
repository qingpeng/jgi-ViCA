import khmer_tst_utils as utils
import subprocess
import os
import shutil

def scriptpath(script):
    return script
    
    
def test_shred_fixed():

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data('example.mito.fasta'), infile)
    
    outfile = infile+'.shred'
    script = scriptpath('shred.py')
    args = ["--input", infile, "--shred", "fixed",     "--out", outfile, 
    "--testing"]
    utils.runscript(script, args, in_dir)
    assert os.path.exists(outfile), outfile
    data = [x.strip() for x in open(outfile)]
    assert len(data) == 6715
    assert data[0].startswith(">gi|408772811|ref|NC_018815.1") == True
    assert data[1] == "CGAGCTACTCCAAGGCAGTCTAATTTATAGGACCCCCCCCGTCTCTGTGGCAAAAGAGTG"
    assert data[84] == "ACATGTTTACTGTAGGAATG"
    assert data[85].startswith(">gi|545906466|ref|NC_022258.1||pos|11965..16965 gi|408772797") == True
    assert data[-1] == "GCCCTTTACGTTACCTGGTC"
    utils.cleanup()
    
    
def test_shred_lognorm():

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data('example.mito.fasta'), infile)
    
    outfile = infile+'.shred'
    script = scriptpath('shred.py')
    args = ["--input", infile, "--shred", "lognorm", "--shape", "1.5","--scale", "1000", 
    "--loc", "3000", "--out", outfile, 
    "--testing"]
    utils.runscript(script, args, in_dir)
    assert os.path.exists(outfile), outfile
    data = [x.strip() for x in open(outfile)]
    assert len(data) == 5706
    assert data[0].startswith(">gi|545906466|ref|NC_022258.1") == True
    assert data[1] == "AGGAAACCTTATTCAATCCGAAGATACCATGAAAATGGGAAATACGGCAAAAGTAAATTA"
    assert data[75] == "GTATTAAAACAAGATTTTTTAACA"
    assert data[76].startswith(">gi|545601253|ref|NC_022435.1|") == True
    assert data[-1] == "GGTTACTCTCGTTGCGGTTATTGCTAAAAGCAAGGGATATTCACTCAGCAGAAAG"
    utils.cleanup()

def test_feature_extraction_metamark():

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data('example.mito.fasta.shreded.subset'), infile)
    
    outfile = infile+'.metamark_vector'
    script = scriptpath('feature_extraction_metamark.py')
    mmp = os.path.abspath("../scripts/gm_parameters/par_11.modified")
    print mmp
    tmp = os.path.abspath("./")
    print tmp
    print in_dir
    args = ["--input", infile, "--outfile", outfile, "--tmp", tmp, "--mmp", 
    mmp, "--taxid", "12345" ]
    utils.runscript(script, args, in_dir)
    assert os.path.exists(outfile), outfile
    print outfile
    data = [x.strip() for x in open(outfile)]
    print len(data)
    assert len(data) == 30
    assert data[1].startswith("12345\tgi|511782593|ref|NC_021399.1") == True
    assert data[0].endswith("0.018679\t0.016415\t0.016415") == True
    assert data[-1].startswith("12345\tgi|511782593|ref|NC_021399.1||pos|295582..300582") == True
    assert data[-1].endswith("0.023325\t0.019296\t0.021310") == True
    utils.cleanup()
    
def test_feature_extraction_kmer():

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data('example.mito.fasta.shreded.subset'), infile)
    
    outfile = infile+'.kmer_vector'
    script = scriptpath('feature_extraction_kmer.py')
    mmp = os.path.abspath("../scripts/gm_parameters/par_11.modified")
    print mmp
    tmp = os.path.abspath("./")
    print tmp
    print in_dir
    args = ["--input", infile, "--outfile", outfile, "--taxid", "12345", "--label", "taxid"]
    utils.runscript(script, args, in_dir)
    assert os.path.exists(outfile), outfile
    print outfile
    data = [x.strip() for x in open(outfile)]
    print len(data)
    assert len(data) == 54
    assert data[1].startswith("12345\tgi|511782593|ref|NC_021399.1") == True
    assert data[0].endswith("0.00300180108065\t0.00640384230538\t0.00380228136882") == True
    assert data[-1].startswith("12345\tgi|511782593|ref|NC_021399.1||pos|304493..309493") == True
    assert data[-1].endswith("0.00300180108065\t0.00500300180108\t0.00200120072043") == True
    utils.cleanup()
    