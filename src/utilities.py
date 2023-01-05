import numpy as np, requests
from bs4 import BeautifulSoup as bs

class custom_iter: # custom iterator class that allows for retrieval of current element w/out advancing
    def __init__(self, iterable):
        self.iterable = iterable
        self.iterator = iter(iterable)
        self.current = None
    def __next__(self):
        try:
            self.current = next(self.iterator)
        except StopIteration:
            self.current = None
        finally:
            return self.current
    def __len__(self):
        return len(self.iterable)

def split_filters(string):
    UVOT_filters = ["B","U","UVW1","UVM2","UVW2","White"]
    name_idxs = custom_iter([string.index(i) for i in UVOT_filters if i in string])
    split_list = [string[name_idxs.current:next(name_idxs)] for i in range(len(UVOT_filters))]
    split_list = [item for item in split_list if len(item)>0]
    return np.unique(split_list).tolist()

def new_since_Fong(dataframe, colname="GRB"):
    indexer = [int(grb[:6]) > 150301 for grb in dataframe[colname]]
    return dataframe[indexer].copy()

def simbad_bibcodes(GRB):
    URL = f"http://simbad.u-strasbg.fr/simbad/sim-id?Ident=GRB%20{GRB}&submit=In+table&output.format=ASCII"
    content = requests.get(URL).text
    entries = content.split("\n\n")
    bibcodes = entries[["Bibcodes" in entry or "References" in entry for entry in entries].index(True)].strip()
    
    return bibcodes.split()[4:]

def literature_references(GRB,titles=True,links=True,GCNs=False):
    title_list = []
    link_list = []
    for bibcode in simbad_bibcodes(GRB):
        if "GCN" in bibcode and GCNs==False:
            continue
        ADS_URL = f"https://ui.adsabs.harvard.edu/abs/{bibcode}/"
        link_list.append(ADS_URL)
        soup = bs(requests.get(ADS_URL).text, features="lxml")
        title = soup.find("title")
        title_list.append(title.text[:-11]) # exclude " - NASA/ADS" from the title
        
    if titles and links:
        return dict(zip(title_list,link_list))
    elif titles:
        return title_list
    elif links:
        return link_list
    else:
        return None
    
def addYear(GRB_df):
    yrs = []
    for i in GRB_df.index:
        yy = int(GRB_df.loc[i,"GRB"][:2])
        if yy>21:
            yrs.append(1900+yy)
        else:
            yrs.append(2000+yy)
    GRB_df["Year"] = yrs