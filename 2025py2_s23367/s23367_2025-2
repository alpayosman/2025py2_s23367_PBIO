from Bio import Entrez, SeqIO
import time,csv,matplotlib.pyplot as plt,pandas as pd
class GB:
 def __init__(s,e,k):Entrez.email,Entrez.api_key=e,k
 def search(s,t):
  with Entrez.esearch(db="nucleotide",term=f"txid{t}[Organism]",usehistory="y")as h:r=Entrez.read(h)
  return r["WebEnv"],r["QueryKey"],int(r["Count"])
 def fetch(s,w,q,c):
  R=[]
  for st in range(0,min(c,2000),500):
   with Entrez.efetch(db="nucleotide",rettype="gb",retmode="text",retstart=st,retmax=500,webenv=w,query_key=q)as h:
    R+=list(SeqIO.parse(h,"genbank"))
   time.sleep(0.4)
  return R
 def filt(s,R,mx,mn):return[r for r in R if mn<=len(r.seq)<=mx]
 def csv(s,R,f):
  with open(f,"w",newline="")as o:
   w=csv.writer(o);w.writerow(["Acc","Len","Desc"])
   [w.writerow([r.id,len(r.seq),r.description]) for r in R]
 def plot(s,R,f):
  d=pd.DataFrame([(r.id,len(r.seq))for r in R],columns=["A","L"]).sort_values("L",ascending=False)
  plt.figure(figsize=(10,5));plt.plot(d["A"],d["L"],marker="o");plt.xticks(rotation=90,fontsize=6)
  plt.tight_layout();plt.savefig(f);plt.close()
if __name__=="__main__":
 e=input("Email: ")
 k=input("API key: ")
 t=input("TaxID: ")
 mn=int(input("Min len: "))
 mx=int(input("Max len: "))
 g=GB(e,k)
 print("Searching...");w,q,c=g.search(t)
 print(f"Found {c}");R=g.fetch(w,q,c)
 F=g.filt(R,mx,mn);print(f"{len(F)} match")
 g.csv(F,f"taxid_{t}.csv")
 g.plot(F,f"taxid_{t}.png")
 print("Done.")
