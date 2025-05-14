import numpy as np
import matplotlib.pyplot as plt
import string
table = str.maketrans('', '', string.ascii_lowercase)


def plot_msa_v2(feature_dict, name,sort_lines=True, dpi=100):
    seq = feature_dict["msa"][0]
    if "asym_id" in feature_dict:
      Ls = [0]
      k = feature_dict["asym_id"][0]
      for i in feature_dict["asym_id"]:
        if i == k: Ls[-1] += 1
        else: Ls.append(1)
        k = i
    else:
      Ls = [len(seq)]
    Ln = np.cumsum([0] + Ls)

    try:
        N = feature_dict["num_alignments"][0]
    except:
        N = feature_dict["num_alignments"]

    msa = feature_dict["msa"][:N]

    gap = msa != '-'
    qid = msa == seq
    #print(gap,qid)
    gapid = np.stack([gap[:,Ln[i]:Ln[i+1]].max(-1) for i in range(len(Ls))],-1)
    lines = []
    Nn = []
    for g in np.unique(gapid, axis=0):
        i = np.where((gapid == g).all(axis=-1))
        qid_ = qid[i]
        gap_ = gap[i]
        seqid = np.stack([qid_[:,Ln[i]:Ln[i+1]].mean(-1) for i in range(len(Ls))],-1).sum(-1) / (g.sum(-1) + 1e-8)
        non_gaps = gap_.astype(float)
        non_gaps[non_gaps == 0] = np.nan
        if sort_lines:
            lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(),None]
        else:
            lines_ = non_gaps[::-1] * seqid[::-1,None]
        Nn.append(len(lines_))
        lines.append(lines_)

    Nn = np.cumsum(np.append(0,Nn))
    lines = np.concatenate(lines,0)
    plt.figure(figsize=(8,5), dpi=dpi)
    plt.title("Sequence coverage"+name)
    plt.imshow(lines,
              interpolation='nearest', aspect='auto',
              cmap="rainbow_r", vmin=0, vmax=1, origin='lower',
              extent=(0, lines.shape[1], 0, lines.shape[0]))
    for i in Ln[1:-1]:
        plt.plot([i,i],[0,lines.shape[0]],color="black")
    for j in Nn[1:-1]:
        plt.plot([0,lines.shape[1]],[j,j],color="black")

    plt.plot((np.isnan(lines) == False).sum(0), color='black')
    plt.xlim(0,lines.shape[1])
    plt.ylim(0,lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")


##############################################################################
    plt.savefig("./plots/"+name+"_coverage.png", bbox_inches='tight')
##############################################################################
    plt.close()



##############################################################################
f=open("compare_out40","w")
##############################################################################


import os

##############################################################################
directory = './dout_sensitive_k0_db0' #gotta be d
##############################################################################

for filename in os.listdir(directory):

##############################################################################
  fname="./dout_sensitive_k0_db0/"+filename #gotta be d
##############################################################################

  file1 = open(fname, 'r')
  Lines = file1.readlines()
  Lines = [i for i in Lines if i != "\x00" and i!="\n"]


  d=[]
  msa=[]
  for i in range(len(Lines)):
    if i%2==1:
      msa.append(Lines[i][:-1])
    elif i!=0:
      d.append(Lines[i][1:])


  if msa:
    m=len(msa[0])


    msa_=[[char for char in string] for string in msa]
    feature_dict={}
    feature_dict['msa']=np.array(msa_)
    feature_dict['asym_id']=[0]*len(msa_[0])
    feature_dict['num_alignments']=[len(msa_)]*len(msa_[0])
##############################################################################
    #A=plot_msa_v2(feature_dict,filename[:-4]+" Diamond")
##############################################################################


##############################################################################
  fname="./mmout1/"+filename
##############################################################################
  file1 = open(fname, 'r')
  Lines = file1.readlines()
  Lines = [i for i in Lines if i != "\x00" and i!="\n"]

  msa=[]
  mm=[]
  for i in range(len(Lines)):
    if i%2==1:
      msa.append(Lines[i][:-1])
    elif i!=0:
      mm.append(Lines[i][1:])

  if msa:
    m=len(msa[0])
    for k in range(len(msa)):
      s=msa[k]
      msa[k]=s.translate(table)

    msa_=[[char for char in string] for string in msa]
    feature_dict={}
    feature_dict['msa']=np.array(msa_)
    feature_dict['asym_id']=[0]*len(msa_[0])
    feature_dict['num_alignments']=[len(msa_)]*len(msa_[0])
##############################################################################
    #B=plot_msa_v2(feature_dict,filename[:-4]+" MMseqs")
##############################################################################


    d=set(d)
    mm=set(mm)
    print(filename[:-4],"\t",len(d),"\t",len(mm),"\t",len(d-mm),"\t",len(mm-d))#,"\t",d-mm,"\t",mm-d,file=f)
f.close()