from Bio import Entrez
import argparse
import pandas as pd
import json

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--email", help = "email")
args = parser.parse_args()

Entrez.email = args.email

infile = "./sa2024_biosample.txt"
indata = pd.read_csv(infile, sep = "\t")
indata.fillna("Not Available", inplace = True)
print(indata)


#collection_date: collection date OR collection_date, collection date
#country: country, geographic location
#iso source: isolation_source, isolation source
#specific host: specific host, host

def searchstringsplit(searchstring, harmname, searchterm):
    split1 = f"harmonized_name=\"{harmname}\""
    split2 = f"display_name=\"{searchterm}\">"
    teststring = searchstring.split(split1)[1].split(split2)[1].split("</Attribute>")[0]
    return teststring
    #outstring = searchstring.split(f"<Attribute attribute_name=\"{searchterm}\">")[1].split("</Attribute>")[0]
    #return outstring    

def getall(searchstring):
    
    #searchstrings = ["collection_date", "geo_loc_name", "isolation_source", "host"]
    #printterms = ["collection date", "geographic location", "isolation source", "host"]

    try:
        collection_date = searchstringsplit(searchstring, "collection_date", "collection date")
    except:
        collection_date = "Not Available"
    
    try:
        geo_loc_name = searchstringsplit(searchstring, "geo_loc_name", "geographic location")
    except:
        geo_loc_name = "Not Available"
    
    try:
        isolation_source = searchstringsplit(searchstring, "isolation_source", "isolation source")
    except:
        isolation_source = "Not Available"
    
    try:
        host = searchstringsplit(searchstring, "host", "host")
    except:
        host = "Not Available"
    
    return collection_date, geo_loc_name, isolation_source, host
    
def main():
    for biosample in indata["Biosample"]:
        handle = Entrez.esearch(db="biosample", term = biosample, retmax = 1)
        record = Entrez.read(handle)
        handle2 = Entrez.esummary(db="biosample", id = record["IdList"][0], retmode = "xml")
        record2 = Entrez.read(handle2)
        searchstring = str(record2["DocumentSummarySet"]["DocumentSummary"][0]["SampleData"])
        
        collection_date, geo_loc_name, isolation_source, host = getall(searchstring)
        indata.loc[indata["Biosample"] == biosample, "Date"] = collection_date
        indata.loc[indata["Biosample"] == biosample, "Geo"] = geo_loc_name
        indata.loc[indata["Biosample"] == biosample, "Source"] = isolation_source
        indata.loc[indata["Biosample"] == biosample, "Host"] = host
        
    indata.to_csv("./sa2024_biosample_meta.txt", sep = "\t", index = False)
main()

    
    
