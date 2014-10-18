import urllib2
from bs4 import BeautifulSoup 
import networkx as nx
import math
from scipy import integrate

def removeNonAscii(s): return "".join(i for i in s if ord(i)<128)

def replace_list (s, replacements):
    return reduce (lambda t, r: t.replace (*r), replacements, s)

def get_page (page_url):
    print page_url
    return BeautifulSoup(urllib2.urlopen(page_url))

def get_name (page):
     try:
         name = page.h1.img['title']
     except:
         name = page.h1.get_text()
     return replace_list (removeNonAscii(name), [('"',''), ("'", "")]) #sometimes there's a double quote in the font name that will destroy your graph

def get_tags (page):
    tags = map (lambda t: replace_list (t.get_text().lower(), [('-', ''), (' ','')]),\
                    page('span', class_ = 'tag'))
    return filter (lambda tag: len(tag)>3, tags)

def name_tags (page):
    try:
        return (get_name(page), get_tags(page))
    except:
        print 'Error extracting informations'
        return False

def get_fonts_with_tags (page_urls):
    return filter( lambda f: f and (len (font_tags (f)) > 0), 
                   map (lambda p: name_tags (get_page (p)), page_urls))


def get_unique_tags (fonts):
    return list(set (reduce (lambda a, b: a+b, map (font_tags, fonts))))

    
def font_name (f): return f[0]

def font_tags (f): return f[1]

def tag_in_font (tag, font):
    return (tag in  font_tags(font))

#this function returns a list of unique tags associated with the fonts they describe
def tag_with_fonts (fonts):
    utags = get_unique_tags (fonts)
    return map (lambda tag: \
                    (tag, map (font_name, filter (lambda font: tag_in_font (tag, font), fonts))),\
                    utags)

def bipartite_network (fonts):
    return (fonts, tag_with_fonts(fonts))

def reverse_bipartite (bip):
    return (bip[1], bip[0])

def tag_name (t): return t[0]

def tag_fonts (t): return t[1]

def name (n): return n[0]
def properties (n): return n[1]


def freq (bip, prop):
    return len (dict (reverse_bipartite (bip)[0])[prop])

def freq_a_b (bip, prop_a, prop_b):
    return len (filter (lambda i: (prop_a in properties (i)) and \
                            (prop_b in properties (i)), bip[0]))

def prob (bip, prop):
    return float(freq (bip, prop))/len (bip[0])

def prob_a_b (bip, prop_a, prop_b):
    return float(freq_a_b(bip, prop_a, prop_b))/len (bip[0])

def pmi (bip, a, b):
    p = prob_a_b(bip, a, b)/(prob (bip, a) * prob (bip, b))
    if p > 1:
        return math.log (p)
    else: return 0

def pmi_noweight (bip, a, b):
    if pmi(bip, a, b) > 0: return 1
    else: return 0

    
### how to use it:
# first get the fonts with tags 
# fonts = get_fonts_with_tags (list_of_links_to_myfonts_pages)

# then build a bipartite container:
# bip = bipartite_network (fonts)
#
# now to find the number of time a tag is associated with some font use freq
#
# freq (bip, name_of_a_tag)
#
# to find the pointwise mutual information of a pair of tags:
# pmi (bip, tag1, tag2)
#
# to find the number of time a font is associated with some font you need to reverse the bipartote network:
# pmi(reverse_bipartite(bip), name_of_font)

def table (bip, weight_fn):
    ps = map (name, bip[1])
    return map (lambda n: map (lambda n1: ((n1 != n) and weight_fn (bip, n, n1)) or 0, ps), ps)

def fonts_table (bip, weight_fn):
    return table (reverse_bipartite(bip), weight_fn)

def tags_table (bip, weight_fn):
    return table (bip, weight_fn)


def build_graph (bip, weight_fn, no_weight=False):
    names = map (lambda n: removeNonAscii(name(n)), reverse_bipartite(bip)[0])
    tb = table(bip, weight_fn)
    g = nx.DiGraph()
    for first in range (len (names)):
        for second in range (len (names)):
            if tb[first][second] != 0:
                if no_weight:
                    g.add_edge(first, second)
                else:
                    g.add_weighted_edges_from([(first, second, tb[first][second])])
    
        try:
            g.node[first]['label'] = names [first]
        except:
            g.add_node(first, label=names[first])
    return g


def fonts_graph (bip, weight_fn,  no_weight=False):
    return build_graph (reverse_bipartite(bip), weight_fn, no_weight)

def tags_graph (bip, weight_fn,  no_weight=False):
    return build_graph (bip, weight_fn, no_weight)


####


def build_graphs (fonts, weight_fn, no_weight=False):
    bip = bipartite_network (fonts)
    return (fonts_graph (bip, weight_fn, no_weight), tags_graph (bip, weight_fn, no_weight))

########



def font_in_tag (font, tag):
    return (font in tag_fonts(tag))


def p_in_tags (p, font):
    if (p in font_tags(font)): return 1
    else: return 0

def tags_in_font (tag1, tag2, font):
    if ((tag1 in font_tags(font)) and (tag2 in font_tags(font))):
        return 1
    else: return 0
    

def a_b_edge (a, b):
    if a == b: return 0
    else:
        return reduce (lambda cumul, tag: cumul + p_in_tags (tag, b),\
                           font_tags (a), 0)

def font_pmi_edge (a, b, fonts):
    if a == b: return 0
    tags = float(len (get_unique_tags (fonts)))
    prob = (a_b_edge(a, b)/tags)/((len (font_tags(a))/tags)*(len (font_tags(b))/tags))
    if prob > 1:
        return math.log(prob)
    else: return 0

def normalized_edge (edge_weight, font):
    return float(edge_weight)/len(font_tags(font))

def table_old (fonts):
    return map (lambda font: \
                    map (lambda font1:\
                             a_b_edge (font, font1), fonts), fonts)

def normalized_table (fonts):
    return map (lambda font: \
                    map (lambda font1: \
                             float(a_b_edge (font, font1))/len(font_tags(font)), \
                             fonts), fonts)

def fonts_pmi_table (fonts):
    return map (lambda font: \
                    map (lambda font1:\
                             font_pmi_edge (font, font1, fonts), fonts), fonts)

def build_fonts_graph (fonts, table_fn):
    names = map (font_name, fonts)
    tb = table_fn(fonts)
    g = nx.DiGraph()
    for first in range (len (names)):
        for second in range (len (names)):
            if tb[first][second] != 0:
                g.add_weighted_edges_from([(first, second, tb[first][second])])
    
        try:
            g.node[first]['label'] = names [first]
        except:
            g.add_node(first, label=names[first])
    return g
                

def tags_edge (tag1, tag2, fonts):
    if tag1==tag2: return 0;
    else:
        return reduce (lambda cumul, font: cumul + tags_in_font (tag1, tag2, font), fonts, 0)
        
def tag_table (utags, fonts):
    return map (lambda tag: \
                    map (lambda tag1:\
                             tags_edge (tag, tag1, fonts), utags), utags)

def prob_a_b_edge (tag1, tag2, fonts):
    return tags_edge(tag1,tag2,fonts)/float (len (fonts))

def prob_tag (tag, fonts):
    return sum (map (lambda f: p_in_tags (tag, f), fonts))/float(len (fonts))

def pmi_tag (tag1, tag2, fonts):
    p = prob_a_b(tag1, tag2, fonts)/(prob_tag (tag1, fonts) * prob_tag (tag2, fonts))
    if p>1:
        return math.log(p)
    else: return 0

def tag_pmi_table (utags, fonts):
    return map (lambda tag: \
                    map (lambda tag1: pmi_tag(tag, tag1, fonts), utags), utags)

def tag_freq (fonts):
    utags = get_unique_tags(fonts)
    freqs = map (lambda tag: \
                     len (filter (lambda f:\
                                      (tag in font_tags(f)), fonts)),\
                     utags)
    return zip (utags, freqs)

def build_tags_graph (fonts, tab_fn):
    utags = get_unique_tags (fonts)
    names = map (lambda n: str (removeNonAscii (n)), utags)
    tb = tab_fn(utags, fonts)
    g = nx.DiGraph()
    for first in range (len (names)):
        for second in range (len (names)):
            if tb[first][second] != 0:
                g.add_weighted_edges_from([(first, second, tb[first][second])])
    
        try:
            g.node[first]['label'] = names [first]
        except:
            g.add_node(first, label=names[first])
    return g
               
def filter_weights (g, ts):
    filtered_graph = nx.Graph()
    for node in g:
        k_n = len (g[node])
        filtered_graph.add_node(node, label=g.node[node]['label'])
        if k_n > 1:
            for nei in g[node]:
                if g[node][nei]['weight'] > ts:
                    filtered_graph.add_edge(node, nei, weight= g[node][nei]['weight'])
    return filtered_graph

def extract_backbone(g, alpha):
  backbone_graph = nx.Graph()
  for node in g:
      k_n = len(g[node])
      if k_n > 1:
          sum_w = sum( g[node][neighbor]['weight'] for neighbor in g[node] )
          for neighbor in g[node]:
              edgeWeight = g[node][neighbor]['weight']
              pij = float(edgeWeight)/sum_w
              f = lambda x: (1-x)**(k_n-2) 
              alpha_ij =  1 - (k_n-1)*integrate.quad(f, 0, pij)[0] 
              if alpha_ij < alpha: 
                  backbone_graph.add_edge( node,neighbor, weight = edgeWeight)
                  backbone_graph.node[node]['label'] = g.node[node]['label']
                  backbone_graph.node[neighbor]['label'] = g.node[neighbor]['label']
      else:
          backbone_graph.add_edge(node, g[node][0], weight = g[node][0]['weight'])
          backbone_graph.node[node]['label'] = g.node[node]['label']
          backbone_graph.node[g[node](0)]['label'] = g.node[g[node](0)]['label']   
    
  return backbone_graph

ft = [('adobe garamond', ['high contrast', 'medium xh', 'stress', 'open counters', 'small eye', 'no spurs', 'serif', 'bracket']),\
          ('bodoni', ['very high contrast', 'small xh', 'no stress', 'geometric', 'ball terminals', 'closed counters', 'medium eye', 'spurs', 'serif', 'straight']),\
          ('baskerville', ['high contrast', 'large xh', 'no stress', 'open counters', 'medium eye' , 'open counters', 'spurs', 'serif' , 'bracket']),\
          ('simoncini garamond', ['high contrast', 'medium xh', 'stress', 'open counters', 'small eye', 'no spurs', 'serif', 'bracket']), \
          ('helvetica', ['no contrast', 'large xh', 'no stress', 'closed counters', 'large eye', 'sans serif', 'horizontal cut', 'mono g', 'a tail']),\
          ('gill sans', ['no contrast', 'medium xh', 'open counters', 'no stress', 'geometric', 'sans serif', 'large eye', 'vertical cut', 'bi g', 'a tail']),\
          ('centaur', ['high contrast', 'small xh', 'open counters', 'stress', 'small eye', 'serif', 'no spurs', 'pen', 'bracket']),\
          ('didot', ['very high contrast', 'small xh', 'closed counters', 'spurs', 'medium eye', 'ball terminals', 'no stress', 'serif', 'straight', 'geometric']),\
          ('jannon', ['high contrast', 'open counters', 'variable stress', 'pen', 'serif', 'straight', 'medium xh', 'small eye', 'no spurs']),\
          ('fournier', ['high contrast', 'open counters', 'medium eye', 'no spurs', 'no stress', 'serif', 'straight']),\
          ('futura', ['no contrast', 'small xh', 'closed counters', 'no stress', 'geometric', 'large eye', 'sans serif', 'vertical cut', 'mono g', 'no a tail']),\
          ('balance', ['reversed', 'no contrast','large eye', 'open counters', 'stress', 'sans serif', 'large xh', 'vertical cut', 'mono g', 'no a tail']),\
          ('clarendon', ['serif', 'contrast', 'closed counters', 'no stress', 'slab', 'bracket', 'serif', 'large eye', 'spurs', 'ball terminals']),\
          ('univers', ['no contrast', 'sans serif' , 'closed counters', 'no stress', 'large eye', 'large xh', 'horizontal cut', 'mono g', 'no a tail']),\
          ('scala sans', ['no contrast', 'sans serif', 'some stress', 'open counters', 'medium eye','medium xh', 'diagonal cut', 'bi g', 'a tail' ]),\
          ('akzidenz grot', ['no contrast', 'sans serif', 'no stress', 'large eye', 'closed counters', 'large xh', 'diagonal cut', 'mono g', 'a tail']),\
          ('franklin got', ['no contrast', 'sans serif', 'no stress',  'large eye', 'closed counters', 'large xh', 'diagonal cut', 'bi g', 'no a tail']),\
          ('frutiger', ['no contrast', 'sans serif', 'no stress',  'large eye', 'open counters', 'medium xh', 'vertical cut', 'mono g', 'no a tail']),\
          ('meta', ['some contrast', 'sans serif', 'stress',  'large eye', 'open counters', 'medium xh', 'diagonal cut', 'bi g', 'a tail']),\
          ('theSans', ['no contrast', 'sans serif', 'stress', 'medium eye', 'open counters', 'large xh', 'diagonal cut', 'no a tail', 'bi g']),\
          ('Romain du Roi', ['very high contrast', 'small xh', 'small eye', 'open counters', 'serif', 'bracket', 'no spurs', 'ball terminals', 'no stress']),\
          ('PMN Cecilia', ['serif', 'no contrast', 'slab', 'straight', 'open counters', 'some stress', 'no spurs', 'large eye', 'large xh']),\
          ('Walbaum linotype', ['serif', 'very high contrast', 'closed counters', 'no stress', 'large eye', 'straight', 'no spurs', 'large xh']),\
          ('Rockwell', ['serif', 'straight', 'no contrast', 'large eye', 'spurs', 'large xh', 'diagonal cut', 'no stress', 'closed counters', 'geometric'])]


ft1 = [('adobe garamond', ['serif', 'high contrast', 'stress', 'open counters', 'small eye', 'brackets', 'angled connection']),\
           ('bodoni', ['serif', 'very high contrast', 'double curve R leg', 'closed counters', 'b with serif', 'medium eye', 'balls', 'spurs', 'continous connection']),\
           ('baskerville' , ['serif', 'high contrast', 'open counters', 'medium eye', 'brackets', 'spurs', 'continous connection']),\
           ('helvetica', ['sans', 'no contrast', 'double curve R leg', 'mono g', 'closed counters', 'large eye', 'spurs', 'horizontal cut', 'continous connection']),\
           ('meta', ['sans', 'very low contrast', 'open counters', 'large eye', 'diagonal cut', 'angled connection']), \
           ('centaur', ['serif', 'high contrast', 'stress', 'open counters', 'small eye', 'brackets', 'diagonal e bar', 'spurs', 'angled connection']),\
           ('fournier' , ['serif', 'high contrast', 'double curve R leg', 'open counters', 'b with serif', 'medium eye', 'angled connection']),\
           ('clarendon', ['serif', 'contrast',  'double curve R leg', 'closed counters', 'large eye', 'balls', 'spurs', 'continous connection']),\
           ('didot', ['serif', 'very high contrast', 'double curve R leg', 'closed counters', 'b with serif', 'large eye', 'balls', 'spurs', 'continous connection']),\
           ('scala sans', ['sans', 'stress', 'very low contrast', 'open counters', 'medium eye', 'diagonal cut', 'angled connection']),\
           ('franklin gothic', ['sans', 'very low contrast', 'closed counters', 'large eye', 'spurs', 'diagonal cut', 'continous connection']),\
           ('grotesque mt', ['sans', 'no contrast', 'mono g', 'closed counters', 'large eye', 'diagonal cut', 'continous connection']),\
           ('din', ['sans', 'no contrast', 'mono g', 'open counters', 'large eye', 'geometric', 'diagonal cut', 'continous connection']),\
           ('frutiger', ['sans', 'no contrast', 'mono g', 'open counters', 'large eye', 'vertical cut', 'continous connection']),\
           ('univers', ['sans', 'no contrast', 'mono g', 'large eye', 'horizontal cut', 'closed counters', 'continous connection']),\
           ('jannon', ['serif', 'high contrast', 'stress', 'open counters', 'small eye', 'angled connection']),\
           ('walbaum', ['serif', 'very high contrast', 'closed counters', 'medium eye', 'double curve R leg', 'continous connection']),\
           ('akzidenz', ['sans', 'no contrast', 'closed counters', 'diagonal cut', 'mono g', 'large eye', 'spurs', 'continous connection']),\
           ('fleischmann', ['serif', 'high contrast', 'double curve R leg', 'open counters', 'medium eye', 'brackets', 'spurs', 'balls', 'angled connection']),\
           ('minion', ['serif', 'high contrast', 'open counters', 'medium eye', 'stress', 'brackets', 'angled connection']),\
           ('caslon', ['serif', 'high contrast', 'open counters', 'small eye', 'brackets', 'spurs', 'angled connection']),\
           ('alright sans', ['sans', 'no contrast', 'open counters', 'large eye', 'diagonal cut', 'angled connection']),\
           ('miller', ['serif', 'high contrast', 'double curve R leg', 'closed counters', 'medium eye', 'balls', 'spurs', 'angled connection']),\
           ('eurostile', ['sans', 'no contrast', 'closed counters', 'large eye', 'horizontal cut', 'continous connection', 'square']),\
           ('swift', ['serif', 'hight contrast', 'open counters', 'stress', 'angled connection', 'medium eye']),\
           ('news gothic', ['sans', 'no contrast', 'closed counters', 'large eye', 'spurs', 'diagonal cut', 'angled connection'])
       ]

def filter_tags (fonts, s):
    return map (lambda f: (font_name(f), \
                        filter (lambda tag:  (s not in tag), font_tags(f))), fonts)

serifs_os = ['http://www.myfonts.com/fonts/mti/centaur/', \
              'http://www.myfonts.com/fonts/adobe/stempel-garamond/', \
              'http://www.myfonts.com/fonts/storm/jannon-pro/', \
              'http://www.myfonts.com/fonts/bitstream/goudy-old-style/', \
              'http://www.myfonts.com/fonts/adobe/palatino/', \
              'http://www.myfonts.com/fonts/linotype/adobe-garamond/', \
              'http://www.myfonts.com/fonts/linotype/sabon/', \
              'http://www.myfonts.com/fonts/linotype/granjon/', \
              'http://www.myfonts.com/fonts/linotype/simoncini-garamond/', \
              'http://www.myfonts.com/fonts/mti/bembo/', \
              'http://www.myfonts.com/fonts/adobe/dante/',\
             'http://www.myfonts.com/fonts/adobe/caslon/',\
                 'http://www.myfonts.com/fonts/itc/galliard/', \
                'http://www.myfonts.com/fonts/typofonderie/le-monde-livre-classic/',\
                 'http://www.myfonts.com/fonts/fontfont/ff-scala/',\
                 'http://www.myfonts.com/fonts/underware/dolly/',\
                 'http://www.myfonts.com/fonts/mti/plantin/',\
                 'http://www.myfonts.com/fonts/adobe/itc-giovanni/',\
                 'http://www.myfonts.com/fonts/adobe/minion/',\
                 'http://www.myfonts.com/fonts/typofonderie/apolline/',\
                 'http://www.myfonts.com/fonts/adobe/garamond-3/',\
                 'http://www.myfonts.com/fonts/mvbfonts/mbv-verdigris-pro/']

serifs_mod = ['http://www.myfonts.com/fonts/linotype/bodoni/', \
                  'http://www.myfonts.com/fonts/urw/firmin-didot/', \
                  'http://www.myfonts.com/fonts/emigre/filosofia-ot/', \
                  'http://www.myfonts.com/fonts/bagraphics/torino-modern/', \
                  'http://www.myfonts.com/fonts/shinn/walburn/', \
                  'http://www.myfonts.com/fonts/mti/walbaum-mt/', \
                  'http://www.myfonts.com/fonts/fontbureau/benton-modern-display/']

trans = ['http://www.myfonts.com/fonts/fontbureau/whitman/',\
         'http://www.myfonts.com/fonts/bitstream/charter/', \
             'http://www.myfonts.com/fonts/emigre/mrs-eaves-ot/', \
             'http://www.myfonts.com/fonts/fontbureau/farnham/', \
             'http://www.myfonts.com/fonts/mti/baskerville-mt/', \
             'http://www.myfonts.com/fonts/adobe/utopia/', \
             'http://www.myfonts.com/fonts/fontbureau/miller/',\
             'http://www.myfonts.com/fonts/adobe/new-caledonia/',\
             'http://www.myfonts.com/fonts/adobe/bulmer/']

slabs = ['http://www.myfonts.com/fonts/typofonderie/le-monde-courrier-std/',\
             'http://www.myfonts.com/fonts/exljbris/museo-slab/', \
             'http://www.myfonts.com/fonts/adobe/clarendon/', \
             'http://www.myfonts.com/fonts/itc/century/']


sans = ['http://www.myfonts.com/fonts/marksimonson/proxima-nova/',\
    'http://www.myfonts.com/fonts/hvdfonts/brandon-grotesque/',\
    'http://www.myfonts.com/fonts/fontbureau/interstate/',\
    'http://www.myfonts.com/fonts/adobe/helvetica-neue/',\
    'http://www.myfonts.com/fonts/bitstream/futura/',\
    'http://www.myfonts.com/fonts/linotype/univers/',\
    'http://www.myfonts.com/fonts/processtype/klavika/',\
    'http://www.myfonts.com/fonts/bitstream/news-gothic/',\
    'http://www.myfonts.com/fonts/berthold/akzidenz-grotesk-bq/',\
    'http://www.myfonts.com/fonts/adobe/gill-sans/',\
    'http://www.myfonts.com/fonts/fontbureau/bureau-grot/',\
    'http://www.myfonts.com/fonts/fontfont/ff-meta/',\
    'http://www.myfonts.com/fonts/fontfont/ff-scala-sans/',\
            'http://www.myfonts.com/fonts/type-together/bree/',\
            'http://www.myfonts.com/fonts/fontbureau/fb-nobel/',\
            'http://www.myfonts.com/fonts/type-together/adelle-sans/',\
            'http://www.myfonts.com/fonts/urw/alternate-gothic/',\
            'http://www.myfonts.com/fonts/linotype/eurostile/',\
            'http://www.myfonts.com/fonts/itc/franklin-gothic/',\
            'http://www.myfonts.com/fonts/mti/grotesque-mt/grotesque-mt/']

links = ['http://www.myfonts.com/fonts/marksimonson/proxima-nova/', 'http://www.myfonts.com/fonts/hvdfonts/brandon-grotesque/', 'http://www.myfonts.com/fonts/fenotype/mercury-script/', 'http://www.myfonts.com/fonts/font-fabric/nexa/', 'http://www.myfonts.com/fonts/sudtipos/hipster-script-pro/', 'http://www.myfonts.com/fonts/exljbris/museo-sans/', 'http://www.myfonts.com/fonts/yellow-design/verb/','http://www.myfonts.com/fonts/berthold/akzidenz-grotesk-be/', 'http://www.myfonts.com/fonts/emily-lime/bombshell-pro/', 'http://www.myfonts.com/fonts/bitstream/futura/', 'http://www.myfonts.com/fonts/yellow-design/veneer/', 'http://www.myfonts.com/fonts/emily-lime/carolyna-pro-black/', 'http://www.myfonts.com/fonts/rene-bieder/rbno2-1/', 'http://www.myfonts.com/fonts/exljbris/museo-slab/', 'http://www.myfonts.com/fonts/cultivated-mind/luella/', 'http://www.myfonts.com/fonts/linotype/avenir/', 'http://www.myfonts.com/fonts/hvdfonts/pluto-sans/', 'http://www.myfonts.com/fonts/hvdfonts/diamonds/', 'http://www.myfonts.com/fonts/linotype/univers/', 'http://www.myfonts.com/fonts/bitstream/swiss-721/', 'http://www.myfonts.com/fonts/paratype/pt-sans-pro/', 'http://www.myfonts.com/fonts/emtype/geogrotesque/', 'http://www.myfonts.com/fonts/hvdfonts/pluto/', 'http://www.myfonts.com/fonts/parachute/pf-din-text-pro/', 'http://www.myfonts.com/fonts/burodunst/novel-sans-rounded-pro/', 'http://www.myfonts.com/fonts/processtype/klavika/', 'http://www.myfonts.com/fonts/font-fabric/intro/', 'http://www.myfonts.com/fonts/linotype/trade-gothic/', 'http://www.myfonts.com/fonts/latinotype/ride-my-bike/', 'http://www.myfonts.com/fonts/corradine/neuron/', 'http://www.myfonts.com/fonts/hubertjocham/narziss/', 'http://www.myfonts.com/fonts/pintassilgo/populaire/', 'http://www.myfonts.com/fonts/fenotype/mishka/', 'http://www.myfonts.com/fonts/bitstream/clarendon/','http://www.myfonts.com/fonts/exljbris/museo/', 'http://www.myfonts.com/fonts/fontbureau/benton-sans/', 'http://www.myfonts.com/fonts/okay-type/alright-sans/', 'http://www.myfonts.com/fonts/parachute/pf-din-text-cond-pro/', 'http://www.myfonts.com/fonts/linotype/avenir-next-pro/', 'http://www.myfonts.com/fonts/correspondence-ink/belluccia/', 'http://www.myfonts.com/fonts/argentina-lian-types/aire/', 'http://www.myfonts.com/fonts/type-together/adelle/', 'http://www.myfonts.com/fonts/emily-lime/carolyna/', 'http://www.myfonts.com/fonts/parachute/pf-centro-sans-pro/', 'http://www.myfonts.com/fonts/linotype/frutiger/', 'http://www.myfonts.com/fonts/font-fabric/code-pro/', 'http://www.myfonts.com/fonts/linotype/neue-helvetica/']

fonts = [('Proxima Nova', [u'clean', u'geometric', u'industrial sans', u'informal', u'sans-serif', u'superfamily', u'workhorse']), ('Brandon Grotesque', [u'1920s', u'1930s', u'art deco', u'bauhaus', u'classic', u'clean', u'display', u'geometric', u'gothic', u'grotesque', u'legible', u'lineal', u'modern', u'monoline', u'retro', u'sans', u'sans-serif', u'succinct', u'thick']), ('Mercury Script', [u'50s', u'alternates', u'brush', u'calligraphy', u'connected', u'delicate', u'elegant', u'fancy', u'hand drawn', u'lettering', u'mishka', u'organic', u'ornaments', u'script', u'slim tony', u'stylish', u'swash', u'upright', u'verner', u'vintage']), ('Nexa', [u'1950s', u'1960s', u'1980s', u'american', u'art', u'avant garde', u'bauhaus', u'circle', u'circled', u'classic', u'clean', u'contemporary', u'cool', u'cursive', u'display', u'elegant', u'extended', u'fancy', u'futura', u'geometric', u'gothic', u'graphic', u'grotesk', u'grotesque', u'headline', u'italics', u'legible', u'logo', u'lubalin', u'modern', u'polish', u'polski', u'sans', u'sans-serif', u'sans serif', u'signage', u'square', u'swiss', u'text', u'trendy', u'true italic', u'wide']), ('Hipster Script Pro', [u'1940s', u'1950s', u'1960s', u'advertising', u'brush', u'casual', u'commercial', u'contemporary', u'cursive', u'decorative', u'fancy', u'fashionable', u'fast', u'free', u'handmade', u'handwriting', u'headlines', u'hipster', u'illustration', u'informal', u'lettering', u'magazine', u'modern', u'retro', u'showcard', u'sign painting', u'speed', u'speedletter', u'swashes', u'tdc', u'wishlist']), ('Museo Sans', [u'canon', u'capital sharp s', u'clean', u'contemporary', u'corporate', u'cufon', u'elegant', u'fashionable', u'free fonts', u'geometric', u'grotesk', u'information', u'international style', u'legible', u'linear', u'linear sans', u'magazine', u'modern', u'modernism', u'modest', u'neo-grotesque', u'neutral', u'news', u'newsletter', u'sans-serif', u'signage', u'simple', u'sturdy', u'technical', u'text', u'versal eszett', u'workhorse']), ('Verb', [u'black', u'bold', u'business', u'clean', u'contemporary', u'corporate', u'display', u'free', u'geometric', u'headline', u'heavy', u'industrial sans', u'italics', u'legible', u'light', u'magazine', u'metro', u'modern', u'oldstyle', u'open', u'poster', u'sans', u'sans-serif', u'sans serif', u'small caps', u'sports', u'sturdy', u'tabular', u'technical', u'text', u'thick', u'thin', u'true italics', u'typographic', u'typography', u'workhorse']), (u'Akzidenz-Grotesk\ufffd BE', [u'gothic', u'grotesk', u'legible', u'modern', u'modernism', u'realist', u'swiss']), ('Bombshell Pro', [u'alternates', u'beautiful', u'branding', u'calligraphic', u'calligraphy', u'connected', u'contrast', u'copperplate', u'cursive', u'decorative', u'display', u'elegant', u'fashion', u'feminine', u'girly', u'handwriting', u'handwritten', u'headline', u'high contrast', u'invitation', u'joined', u'legible', u'ligatures', u'logo', u'love', u'nib', u'original', u'passionate', u'penmanship', u'quirky', u'romantic', u'rough', u'script', u'sexy', u'signature', u'swashes', u'tattoo', u'unique', u'wedding']), (u'Futura\ufffd', [u'1920s', u'1980s', u'abs-cbn', u'advertising', u'animax', u'bauhaus', u'classic', u'geometric', u'julie-sans', u'legible', u'lineal', u'magazines', u'marketing', u'round', u'sans-serif']), ('Veneer', [u'1900s', u'all caps', u'american', u'grotesque', u'grunge', u'industrial', u'letterpress', u'retro', u'sans serif', u'stamp', u'vintage', u'weathered', u'wood', u'wood type', u'worn']), ('Carolyna Pro Black', [u'alternates', u'bold', u'calligraphic', u'calligraphy', u'chase your dreams', u'connected', u'contrast', u'copperplate', u'curly', u'cursive', u'decorative', u'fancy', u'fashion', u'fat', u'feminine', u'formal', u'friendly', u'fully-connected', u'girly', u'hand', u'hand-drawn', u'handwritten', u'high contrast', u'invitation', u'joined', u'legible', u'lettering', u'ligatures', u'nib', u'original', u'pen', u'pretty', u'quirky', u'swash', u'swashes', u'wedding', u'whimsical']), ('RBNo2.1', [u'condensed', u'geometric', u'modern', u'remember_sans', u'sans', u'squarish', u'technical']), ('Museo Slab', [u'capital sharp s', u'cufon', u'display', u'friendy', u'geometric', u'headline', u'legible', u'magazine', u'robust', u'slab', u'slab serif', u'text', u'uec-script', u'versal eszett', u'wishlist']), ('LUELLA', [u'bold', u'book', u'borders', u'candy', u'casual', u'catchwords', u'children', u'cute', u'deco', u'decorative', u'dingbats', u'doodles', u'elegant', u'feminine', u'flowers', u'frames', u'fun', u'funny', u'greeting card', u'hand', u'hand-drawn', u'hand-made', u'handmade', u'handwriting', u'heart', u'illustrations', u'invitation', u'kids', u'labels', u'leaves', u'magazine', u'menu', u'modern', u'opentype', u'ornaments', u'pen', u'retro', u'tall', u'uneven', u'vintage', u'wedding']), ('Avenir', [u'at&t', u'clean', u'elegant', u'geometric', u'legible', u'magazine', u'sans-serif', u'swiss', u'xsf']), ('Pluto Sans', [u'alternates', u'arrows', u'classic', u'clean', u'conceptual', u'condensed', u'cute', u'elegant', u'extended', u'family', u'fashionable', u'fractions', u'friendly', u'futuristic', u'geometric', u'german', u'grotesk', u'high x-height', u'hvd', u'legible', u'linear', u'magazine', u'modern', u'monoline', u'narrow', u'neutral', u'readable', u'sans', u'sans-serif', u'sci-fi', u'signage', u'soft', u'static', u'superfamily', u'super family', u'tabular', u'technical', u'technology', u'text', u'web', u'webfont', u'workhorse']), ('Diamonds', [u'all caps', u'alternates', u'arrows', u'avant garde', u'black', u'circle', u'clean', u'cool', u'cosmetics', u'display', u'elegant', u'experimental', u'fashionable', u'feminine', u'futura', u'geometric', u'grotesk', u'headline', u'hipster', u'legible', u'ligatures', u'light', u'linear', u'logo', u'luxury', u'magazine', u'modern', u'monoline', u'ornaments', u'retro', u'sans', u'sans-serif', u'sans serif', u'simple', u'technical', u'thin', u'web', u'webfont']), ('Univers', [u'business text', u'grotesk', u'information', u'legible', u'magazine', u'news', u'plain', u'sans-serif', u'signage', u'swiss', u'tech pubs', u'traffic', u'transport', u'univers-like']), ('Swiss 721', [u'1950s', u'1980s', u'business text', u'corporate', u'german', u'grotesk', u'helvetica', u'information', u'legible', u'linear', u'linear sans', u'magazine', u'modern', u'modest', u'newsletter', u'sans-serif', u'signage', u'static', u'swiss', u'technical', u'tech pubs', u'xsf']), ('PT Sans Pro', [u'banner', u'black', u'bold', u'business', u'business text', u'capital sharp s', u'caption', u'clean', u'clear', u'commercial', u'condensed', u'corporate', u'cyrillic', u'easy', u'editorial', u'family', u'highway', u'humanist', u'information', u'legible', u'ligatures', u'light', u'modern', u'modest', u'neutral', u'news', u'newspaper', u'optical sizes', u'plain', u'professional', u'regular', u'sans-serif', u'sans serif', u'signage', u'small caps', u'street', u'subway', u'superfamily', u'tech pubs', u'text', u'traffic', u'transport', u'true italics', u'turkish', u'versal eszett', u'workhorse']), ('Geogrotesque', [u'clean', u'din', u'editorial', u'fresh', u'grotesk', u'grotesque', u'magazine', u'neutral', u'publications', u'rounded', u'sans-serif', u'sansserif', u'signage', u'square', u'squarish', u'stock-sans', u'stockmarr', u'techno', u'text']), ('Pluto', [u'cafe', u'cool', u'cupcake', u'curly', u'cutesy', u'decorative', u'elegant', u'emblem', u'fashionable', u'feminine', u'feminine headline', u'fluid', u'french', u'friendly', u'geometric', u'girly', u'graceful', u'headline', u'high x-height', u'legible', u'like_it', u'retro', u'rounded serif', u'sans-serif', u'scripty sans', u'soft', u'text', u'upright', u'upright italic', u'upright script']), (u'PF Din Text Pro\ufffd', [u'corporate', u'cyrillic', u'din', u'grotesk', u'lacks u+00a0', u'legible', u'magazine', u'sans-serif', u'signage', u'superfamily', u'true italics']), ('Novel Sans Rounded Pro', [u'branding', u'business', u'contemporary', u'corporate', u'elegant', u'food', u'humanist', u'information', u'lacks u+00a0', u'linear', u'magazin', u'magazine', u'modern', u'office', u'pachaging', u'round', u'rounded', u'sans', u'soft', u'superfamily', u'super family', u'text', u'web', u'workhorse']), ('Klavika', [u'contemporary', u'facebook', u'geometric', u'modern', u'sans', u'sans-serif', u'splayed m', u'square', u'squarish', u'straight-leg r', u'text', u'true italics', u'web 2.0']), ('Intro', [u'all caps', u'alternates', u'black', u'bulgarian', u'capitals', u'caps', u'classic', u'clean', u'contemporary', u'cool', u'corporate', u'cufon', u'czech', u'elegant', u'fashionable', u'geometric', u'geometry', u'german', u'graceful', u'grotesk', u'grotesque', u'humanist', u'informal', u'inline', u'italic', u'legible', u'light', u'linear', u'linear sans', u'logo', u'magazine', u'modern', u'newsletter', u'polish', u'polski', u'retro', u'romanian', u'russian', u'sans', u'sans-serif', u'sans serif', u'signage', u'spanish', u'square', u'squarish', u'sturdy', u'super family', u'swash', u'technical', u'techno', u'text', u'thin', u'true italic', u'true italics', u'turkish', u'workhorse', u'\u010de\u0161tina', u'\u0440\u0443\u0441\u0441\u043a\u0438\u0439', u'\u0431\u044a\u043b\u0433\u0430\u0440\u0441\u043a\u0438']), ('Trade Gothic', [u'commercial', u'grotesk', u'legible', u'magazine', u'xsf']), ('Ride my Bike', [u'alternates', u'arrow', u'artisan', u'bike', u'bikes', u'brands', u'chilean typefaces', u'condensed', u'coto mendoza', u'cupcake', u'curly', u'decorative', u'delicated', u'dingbats', u'display', u'editorial', u'fashion', u'frames', u'funny', u'guisela mendoza', u'handmade', u'handwriting', u'headlines', u'heart', u'hipster', u'indie', u'indie style', u'informal', u'italic', u'lettering', u'ligatures', u'light', u'logos', u'magazine', u'ornaments', u'retro', u'ridemybike', u'sans serif', u'script', u'street style', u'swash', u'swashes', u'tags', u'workhorse']), ('Neuron', [u'advertising', u'capital sharp s', u'clean', u'contemporary', u'corporate', u'cyrillic', u'friendly', u'headline', u'humanist', u'identity', u'industrial', u'informal', u'label', u'legible', u'logo', u'magazine', u'masculine', u'metro', u'minimal', u'modern', u'neutral', u'organic', u'packaging', u'readable', u'rounded', u'russian', u'sans-serif', u'signage', u'sleek', u'smooth', u'soft', u'south america', u'square', u'tag', u'techno', u'technology', u'true italics', u'urban', u'versal eszett', u'versatile', u'web']), ('Narziss', [u'ball terminals', u'decorative', u'display', u'elegant', u'exclusive', u'fashion', u'hairline', u'hairlines', u'high-contrast', u'serif', u'swash', u'swirls', u'tendrils', u'thick and thin', u'uec-logotype']), ('Populaire', [u'60s', u'1960s', u'art', u'atelier populaire', u'book cover', u'cartoon', u'casual', u'condensed', u'cute', u'decorative', u'doodles', u'friendly', u'funny', u'gothic', u'hand-drawn', u'handcrafted', u'handmade', u'handwriting', u'handwritten', u'illustration', u'informal', u'irregular', u'lettering', u'lively', u'narrow', u'nonprofit', u'organic', u'ornaments', u'party', u'picture', u'random', u'revolt', u'revolution', u'sans-serif', u'script', u'spontaneous', u'ultra-narrow', u'uneven', u'want']), ('Mishka', [u'advertising', u'alternates', u'barber', u'brush', u'connected', u'delicate', u'elegant', u'fancy', u'flourish', u'fresh', u'hand drawn', u'lachrymal terminals', u'lettering', u'mercury script', u'organic', u'ornaments', u'pepita', u'sailor', u'script', u'soft', u'stylish', u'swash', u'tattoo font', u'upright script', u'verner']), ('Clarendon', [u'1800s', u'1960s', u'egyptian', u'legible', u'realist', u'serif', u'slab serif', u'sony']), ('Museo', [u'balanced', u'capital sharp s', u'clean', u'cufon', u'fashionable', u'geometric', u'headline', u'hybrid', u'legible', u'modern', u'monoline', u'semi-serif', u'simple', u'slab serif', u'some-serifs', u'tips and tools', u'unique']), (u'Benton Sans\ufffd', [u'american', u'clean', u'compressed', u'condensed', u'corporate', u'extensive', u'extra', u'favorites', u'gothic', u'grotesque', u'huisstijl', u'iconic', u'leestekst', u'light', u'lineaar', u'magazine', u'mega', u'narrow', u'news', u'newspaper', u'sans-serif', u'standard', u'succinct', u'super', u'thin', u'ultra', u'wide', u'workhorse']), ('Alright Sans', [u'capital sharp s', u'clean', u'conservative', u'contemporary', u'cool', u'corporate', u'editorial', u'exhibitions', u'fashion', u'fashionable', u'friendly', u'grotesk', u'grotesque', u'heavy-hitters', u'humanist', u'industrial sans', u'information', u'legible', u'magazine', u'minimal', u'modern', u'news', u'newspaper', u'sans-serif', u'sans serif', u'signage', u'small caps', u'text', u'thin', u'versal eszett', u'versatile', u'warm', u'workhorse']), (u'PF Din Text Condensed Pro\ufffd', [u'condensed', u'corporate', u'cyrillic', u'din', u'elegant', u'greek', u'grotesk', u'lacks u+00a0', u'legible', u'magazine', u'narrow', u'sans-serif', u'signage', u'superfamily', u'true italics']), ('Belluccia', [u'alternates', u'antiqued', u'artistic', u'brush', u'calligraphic', u'calligraphy', u'classic', u'connected', u'cursive', u'decorative', u'distressed', u'elegant', u'fancy', u'fashion', u'feminine', u'flourishes', u'formal', u'friendly', u'fun', u'girly', u'hand-drawn', u'handwriting', u'handwritten', u'italic', u'ke-collar', u'ke-deluxe-moshup', u'lettering', u'love', u'organic', u'ornaments', u'rough', u'scrapbook', u'shabby chic', u'spencerian', u'swash', u'swash caps', u'tuscan', u'unique', u'unusual', u'vintage', u'wedding', u'whimsical']), ('Aire', [u'1800s', u'alternates', u'beautiful', u'birthday', u'bodoni', u'book', u'book-cover', u'calligraphy', u'card', u'catwalk', u'christmas', u'classic', u'coffee', u'coquette', u'cosmetic', u'cyrillc', u'dazzling', u'dear', u'delicate', u'didone', u'didot', u'display', u'elegant', u'family', u'fancy', u'fashion', u'feminin', u'fleurons', u'flourish', u'flower', u'formal', u'honey', u'invitation', u'italic', u'lian', u'love', u'marriage', u'model', u'modern', u'ornament', u'roman', u'spring', u'stationery', u'sweet', u'thin', u'upright', u'wedding']), ('Adelle', [u'book', u'contemporary', u'eastern european', u'editorial', u'egyptian', u'elegant', u'fashionable', u'fine crafted', u'fluid', u'graceful', u'humanist', u'language support', u'legible', u'low-contrast', u'magazine', u'newspaper', u'oldstyle numerals', u'screen', u'serif', u'slab', u'slab serif', u'smallcaps', u'text', u'title', u'versatile', u'vivid', u'webfont', u'workhorse']), ('Carolyna', [u'alternates', u'calligraphic', u'calligraphy', u'classic', u'connected', u'copperplate', u'curly', u'cursive', u'curvy', u'elegant', u'fancy', u'feminine', u'formal', u'friendly', u'fully-connected', u'girly', u'hand-drawn', u'handwriting', u'handwritten', u'invitation', u'joined', u'joined-up handwriting', u'legible', u'lettering', u'ligatures', u'nib', u'original', u'pretty', u'quirky', u'stylish', u'swash', u'swashes', u'wedding', u'whimsical']), ('PF Centro Sans Pro', [u'clear', u'corporate', u'cyrillic', u'greek', u'humanist', u'legible', u'magazine', u'newspaper', u'sans-serif', u'superfamily', u'text', u'true italics', u'wishlist', u'workhorse']), ('Frutiger', [u'1970s', u'1980s', u'information', u'legible', u'magazine', u'sans-serif', u'signage', u'traffic', u'transport', u'xsf']), ('Code Pro', [u'1960s', u'allcaps', u'all caps', u'capital sharp s', u'caps', u'circle', u'classic', u'clean', u'cool', u'cyrilic', u'disco', u'fonts', u'free', u'futura', u'geometric', u'monoline', u'newspaper', u'party', u'pop', u'sans-serif', u'sans serif', u'uec-sans-geometric', u'versal eszett']), ('Neue Helvetica', [u'1980s', u'business text', u'german', u'grotesk', u'linear', u'modern', u'modest', u'neutral', u'sans-serif', u'static', u'swiss', u'technical', u'tech pubs', u'xsf'])]
