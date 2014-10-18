import networkx as nx
import math
from scipy import integrate

def removeNonAscii(s): return "".join(i for i in s if ord(i)<128)

def replace_list (s, replacements):
    return reduce (lambda t, r: t.replace (*r), replacements, s)

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

def filter_tags (fonts, s):
    return map (lambda f: (font_name(f), \
                        filter (lambda tag:  (s not in tag), font_tags(f))), fonts)



typefaces = [('adobe garamond', ['serif', 'high contrast', 'stress', 'open counters', 'small eye', 'brackets', 'angled connection']),\
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


f, t = build_graphs (typefaces, pmi)
nx.write_gml(t, 'tags.gml')
nx.write_gml(f, 'typefaces.gml')



