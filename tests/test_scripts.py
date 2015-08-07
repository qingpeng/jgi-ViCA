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
    args = ["--input", infile, "--shred", "lognorm", "--shape", "1.5","--scale", "1000", "--loc", 
    "3000", "--out", outfile, 
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


