import requests
import sys
import argparse
try:
    from HTMLParser import HTMLParser
except ImportError:
    from html.parser import HTMLParser
from Bio import SeqIO

class URLParser(HTMLParser):
    urlres = ''
    def handle_starttag(self, tag, attrs):
        if tag == 'meta':
            for i, attr in enumerate(attrs):
                if i == 1:
                    self.urlres=str(attr[1]).lstrip("2;URL=")

class BlastParser(HTMLParser):
    ispre= False
    result = ""
    def handle_starttag(self, tag, attrs):
        if tag == "pre":
            self.ispre=True
    def handle_endtag(self, tag):
        if tag == "pre":
            self.ispre=False
    def handle_data(self, data):
        if self.ispre:
            self.result += data + "\n"
def filter_blast(contig, blast_result, identity, evalue, length, alignment):
    blasts = []
    column_name = ["Query", "Subject_ID", "%identity", "alignment_length", "mismatches", \
                   "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit_score"]
    for line in blast_result.rstrip("\n").split("\n"):
        if line.startswith('#'):
            continue
        blast = line.split("\t")
        if len(blast) != 12:
            sys.stderr.write(line+"\n")
        else:
            blast[0]=contig
            resultBlast = {k:v for k,v in zip(column_name, blast)}
            res_evalue = float(resultBlast.get("evalue"))
            res_identity = float(resultBlast.get("%identity"))
            res_alignment = int(resultBlast.get("q.end")) - int(resultBlast.get("q.start"))
            if res_evalue <= evalue and res_identity >= identity and \
               res_alignment >= (alignment*length)/100:
               blasts.append(resultBlast)
    return blasts


def main():
	sample_id=sys.argv[1:][0]
	with open('./%s_ISfinder.txt'%sample_id,'w') as f:
		with open('./%s.fasta'%sample_id,'r') as g:
			for rec in SeqIO.parse(g,'fasta'):

				r = requests.post('https://www-is.biotoul.fr/blast/ncbiIS.php',\
				#r = requests.post('https://www-is.biotoul.fr/blast.php',\
				data = {'title' : 'test', 'seq' : '%s'%rec.seq, 'seqfile' \
										  : '', 'database' : 'ISfindernt', 'prog' : 'blastn', 'blast' \
										  : 'ok', 'alignment' : '7', 'wordsize' : '11','expect' : '10.0' \
										  , 'gapcost' : '5 2'}, verify=False)
				parser = URLParser()
				parser.feed(r.text)
				o = requests.get('https://www-is.biotoul.fr/blast/' + parser.urlres, verify=False)
				parser = BlastParser()
				parser.feed(o.text)
				result = filter_blast(rec.id,parser.result,90,1e-5,30,0.0000001)
				if len(result)>0:
					for i in result:
						f.write('\t'.join(list(i.values()))+'\n')
						
						
if __name__ == '__main__':
    requests.packages.urllib3.disable_warnings()
    main()
