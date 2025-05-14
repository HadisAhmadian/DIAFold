from Bio import SeqIO

def read_fasta_file(file_path):
    print("in read")
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    print("out read")

    return sequences


def div(fname):
  import numpy as np
  file1 = open(fname, 'r')
  Lines = file1.readlines()
  Lines=[l.split() for l in Lines]

  d=dict()

  for l in Lines:
    if l[0] not in d:
      d[l[0]]=[]  
    d[l[0]].append(l[1:])
  
  return d


def cigar_to_alignment_fold(cigar_string, sequence, query) -> str:
    alignment_s = ""
    alignment_q = ""
    index = 0
    s_curr=0
    q_curr=0
    for char in cigar_string:
        if char.isdigit():
            index = index * 10 + int(char)
        else:
            if char == "M":
                alignment_s += sequence[s_curr:s_curr+index]
                alignment_q += query[q_curr:q_curr+index]
                s_curr=s_curr+index
                q_curr=q_curr+index
                index=0
            elif char == "D":
                #alignment_s += sequence[s_curr:s_curr+index]
                #alignment_q += "-"*index
                s_curr=s_curr+index
                index=0
            elif char == "I":
                alignment_q += query[q_curr:q_curr+index]
                alignment_s += "-"*index
                q_curr=q_curr+index
                index=0
            else:
                raise ValueError(f"Invalid character in CIGAR string: {char}")
    if q_curr<len(query):
      alignment_q+=query[q_curr:]
    if len(alignment_s)<len(alignment_q):
      alignment_s+="-"*(len(alignment_q)-len(alignment_s))

    return alignment_q,alignment_s

def match_to_msa_fold(q_main,Lines,outfname):
  import numpy as np

  f=True
  alignments=[]
  for l in Lines:
    s=l[0]
    i=(int(l[2])-1)
    q_tmp=q_main[0:i]
    s_tmp="-"*len(q_tmp)
    q=q_main[i:]
    if q[-1]=="\n":
      q=q[:-1]
    i=(int(l[3])-1)
    s=s[i:]
    a,b=cigar_to_alignment_fold(l[1],s,q)
    if f:
      alignments.append(q_tmp+a)
      f=False
    alignments.append(s_tmp+b)

##########################################################################
  f = open("./dout_All/"+outfname+".a3m", "w")
########################################################################
  i=-1
  for m in alignments:
    if m[-1]!="\n":
      m+="\n"
    if i==-1:
      f.write(">"+outfname+"\n"+m)
      i=i+1
      continue

    f.write(">"+Lines[i][4]+"\n"+m)
    i=i+1

  f.close()

##############################################################################
queries=read_fasta_file("./DBs/DB.fasta")
##############################################################################

##############################################################################
match_per_query=div('./match_sensitive_k0_c1_All')
##############################################################################


for q, m in match_per_query.items():


##############################################################################
  o_f_n=q.split("|")[1]
##############################################################################


  match_to_msa_fold(queries[q],m,o_f_n)