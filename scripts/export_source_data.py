"""Export lossless plotting data to Excel-readable TSV, and check source consistency."""
from pathlib import Path
import ast, csv, json, hashlib
from Bio import SeqIO, Phylo
from Bio.SeqUtils import gc_fraction

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'source_data'
OUT.mkdir(exist_ok=True)

def definitions(path, names):
    tree=ast.parse(path.read_text())
    nodes=[]
    for node in tree.body:
        if isinstance(node,ast.FunctionDef) and node.name in names:
            nodes.append(node)
        elif isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id in names for t in node.targets):
            nodes.append(node)
    namespace={}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),namespace)
    return namespace

def write(name,header,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(header); w.writerows(rows)

e=definitions(ROOT/'scripts/plot_figure_E.py',{'LABELS','CANDIDATES','REFERENCES'})
tree=Phylo.read(ROOT/'data/figure_E/TerL_final.treefile','newick')
tips=[x.name for x in tree.get_terminals()]
assert len(tips)==19 and len(set(tips))==19
seqs={}
for name in ['TerL_all','TerL_aligned','TerL_trimmed']:
    records=list(SeqIO.parse(ROOT/f'data/figure_E/{name}.faa','fasta'))
    assert len(records)==19 and {r.id for r in records}==set(tips),name
    seqs[name]={r.id:str(r.seq) for r in records}
assert {len(s) for s in seqs['TerL_trimmed'].values()}=={552}
write('Figure_E_taxa_sequences.tsv',['tip_order','sequence_id','display_label','group','raw_amino_acids','aligned_amino_acids','trimmed_amino_acids'],[
    [i+1,id,e['LABELS'].get(id,id),'this_study' if id in e['CANDIDATES'] else 'published_reconstructed_reference' if id in e['REFERENCES'] else 'public_homologue',seqs['TerL_all'][id],seqs['TerL_aligned'][id],seqs['TerL_trimmed'][id]] for i,id in enumerate(tips)])
nodes=list(tree.find_clades(order='preorder')); ids={id(n):f'n{i:02d}' for i,n in enumerate(nodes)}
rows=[]
def visit(n,parent='',order=0):
    support=n.name if not n.is_terminal() else ''
    sh,uf=(support.split('/') if support and '/' in support else ('',''))
    rows.append([ids[id(n)],parent,order,n.name if n.is_terminal() else '',n.branch_length if n.branch_length is not None else '',sh,uf,';'.join(t.name for t in n.get_terminals())])
    for i,c in enumerate(n.clades): visit(c,ids[id(n)],i+1)
visit(tree.root)
write('Figure_E_nodes.tsv',['node_id','parent_id','child_order','tip_id','branch_length','SH_aLRT','UFBoot','descendant_tip_ids'],rows)
write('Figure_E_newick.tsv',['format','tree'],[['Newick',(ROOT/'data/figure_E/TerL_final.treefile').read_text().strip()]])

f=definitions(ROOT/'scripts/plot_figure_F.py',{'classify_product','get_short_label','CATEGORY_COLORS'})
rec=SeqIO.read(ROOT/'data/figure_F/Gv142-2_P1.gbk','genbank')
seq=str(rec.seq).upper(); length=len(seq); gc=gc_fraction(rec.seq)
assert length==37522 and abs(gc*100-43.65971963114973)<1e-9
cds=[]; labelled=False
for feature in rec.features:
    if feature.type!='CDS': continue
    product=feature.qualifiers.get('product',['hypothetical protein'])[0]
    category=f['classify_product'](product); label=f['get_short_label'](product)
    if label=='Integrase':
        if labelled: label=None
        labelled=True
    cds.append([feature.qualifiers.get('locus_tag',[''])[0],str(feature.location),int(feature.location.start)+1,int(feature.location.end),feature.location.strand,product,category,f['CATEGORY_COLORS'][category],label or ''])
assert len(cds)==71
write('Figure_F_CDS.tsv',['locus_tag','biopython_location_0based_halfopen','start_1based_inclusive','end_1based_inclusive','strand','product','plot_category','color_hex','plot_label'],cds)
windows=[]
for start in range(0,length,100):
    end=min(start+500,length); frag=seq[start:end]
    if len(frag)<50: continue
    g=frag.count('G'); c=frag.count('C'); local=(g+c)/len(frag)
    windows.append([start,end,(start+end)/2,len(frag),g,c,local,gc,local-gc,(g-c)/(g+c) if g+c else 0])
assert len(windows)==375
write('Figure_F_GC_tracks.tsv',['start_0based','end_exclusive','midpoint_0based','window_length_bp','G_count','C_count','GC_fraction','genome_GC_fraction','GC_deviation','GC_skew'],windows)
write('Figure_F_genome.tsv',['sequence_id','length_bp','GC_fraction','genbank_topology','plot_geometry'],[[rec.id,length,gc,rec.annotations.get('topology','unknown'),'circular display; not evidence of circular molecule']])
write('Figure_F_sequence_chunks.tsv',['sequence_id','start_1based','end_1based','nucleotides'],[[rec.id,s+1,min(s+10000,length),seq[s:s+10000]] for s in range(0,length,10000)])
SeqIO.write(rec,ROOT/'data/figure_F/Gv142-2_P1.fasta','fasta')
result={'taxa':19,'trimmed_alignment_columns':552,'genome_length_bp':length,'genome_gc_percent':gc*100,'CDS':len(cds),'GC_windows':len(windows),'validation':'PASS: data checks only; tree inference is a separate workflow'}
(OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
