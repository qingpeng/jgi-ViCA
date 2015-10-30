import sys

import re
def multiwordReplace(text, wordDic):
    """
    take a text and replace words that match a key in a dictionary with
    the associated value, return the changed text
    """
    rc = re.compile('|'.join(map(re.escape, wordDic)))
    def translate(match):
        return wordDic[match.group(0)]
    return rc.sub(translate, text)

wordDic = {
"Archaea" : "non-virus",
"Bacteria" : "non-virus",
"ssRNAPhage" : "virus",
"ssRNAVirus" : "virus",
"dsRNAPhage" : "virus",
"dsRNAVirus" : "virus",
"dsDNAPhage" : "virus",
"dsDNAVirus" : "virus",
"ssDNAPhage" : "virus",
"ssDNAVirus" : "virus",
"Retroviruses" : "virus",
"Eukaryota" : "non-virus",
"Mitochondrion" : "non-virus",
"Chloroplast" : "non-virus"}

file_in_fh = open(sys.argv[1],'r')
file_out_fh = open(sys.argv[2],'w')

for line in file_in_fh:
    line2 = multiwordReplace(line, wordDic)
    file_out_fh.write(line2)
    
file_in_fh.close()
file_out_fh.close()

