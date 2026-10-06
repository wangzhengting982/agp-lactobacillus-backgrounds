"""Reconstruct portable analysis copies with path-only substitutions.

Original scripts in data/archive remain untouched. This adapter is retained so
every portability change can be inspected and regenerated.
"""
from pathlib import Path
import hashlib, json

ROOT = Path(__file__).resolve().parents[2]
ARCH = ROOT / 'data/archive'
SRC = Path(__file__).resolve().parent
CHANGES = []

def emit(origin, name, replacements):
    text = origin.read_text(encoding='utf-8-sig')
    before = text
    for old, new in replacements:
        if old not in text:
            raise RuntimeError(f'Expected source fragment absent: {origin.name}: {old}')
        text = text.replace(old, new)
    (SRC/name).write_text(text, encoding='utf-8')
    CHANGES.append(dict(source=str(origin.relative_to(ROOT)), portable=name,
                        source_sha256=hashlib.sha256(origin.read_bytes()).hexdigest(),
                        portable_sha256=hashlib.sha256(text.encode()).hexdigest(),
                        replacements=[dict(old=o,new=n) for o,n in replacements]))

PATHS = """ROOT=Path(__file__).resolve().parents[2]
ARCH=ROOT/'data/archive'
RUN=Path(os.environ.get('AGP_DEEPENING_OUTPUT',ROOT/'runs/deepening'))
PREP=Path(os.environ.get('AGP_PREPARATION_ROOT',ROOT/'runs/preparation/project'))
FROZEN=ARCH/'02_统计输入/冻结统计输入'
A=Path(os.environ.get('AGP_ASSOC_ROOT',PREP/'03_分析与结果/association_models'))
D=RUN/'diet_preprocessing';D.mkdir(parents=True,exist_ok=True)
RAW=ARCH/'02_统计输入/六份源数据'
R=RUN;O=RUN/'results';O.mkdir(parents=True,exist_ok=True)
"""
runtime = "sys.path[:0]=['E:/new article/26_核心关联复核与投稿修订/_runtime','E:/new article/37_审稿意见实质修订与复审/修订记录/_runtime']"
oldpaths = "R=Path(__file__).resolve().parents[1];O=R/'完整分析结果';O.mkdir(exist_ok=True)\nA=Path('E:/new article/37_审稿意见实质修订与复审/新增分析复现归档/输入/03_分析与结果/association_models')\nD=Path('E:/new article/AGP全部数据深化评估');RAW=Path('E:/AAGP科学通报/01_原始输入')"
deep = ARCH/'03_饮食与菌群深化/复现代码'
emit(deep/'common.py','common.py',[(runtime,'# Dependencies are supplied by the project environment.'),(oldpaths,PATHS)])
for p in sorted(deep.glob('0*.py')):
    subs=[]
    if p.name=='04_within_person.py':
        subs=[("Path('E:/new article/29_全部增补分析与稿件修订/审计/纵向分析样本_内部索引.tsv')","FROZEN/'纵向分析样本_内部索引.tsv'"),("Path('E:/new article/05_输入数据/派生输入/lactobacillales_genus_counts_all_samples.tsv')","PREP/'05_输入数据/派生输入/lactobacillales_genus_counts_all_samples.tsv'")]
    elif p.name=='05_asv.py':
        subs=[("P=A.parent/'taxonomy_validation'","P=ARCH/'02_统计输入/冻结统计输入'")]
    elif p.name in ['06_core_robustness.py','08_BMI_strata_resolution.py']:
        subs=[("Path('E:/AAGP科学通报/02_模型输入/analysis_matrix_stool_base.tsv')","RAW/'analysis_matrix_stool_base.tsv'")]
    emit(p,p.name,subs)

diet=ARCH/'02_统计输入/饮食预处理'
emit(diet/'01_inventory.py','00a_inventory.py',[
    (runtime,'# Dependencies are supplied by the project environment.'),
    ("R=Path('E:/AAGP科学通报');O=Path(__file__).parent\nP=Path('E:/new article/37_审稿意见实质修订与复审/新增分析复现归档/输入/03_分析与结果/association_models')",PATHS+"\nP=A;O=D"),
    ("files=list((R/'01_原始输入').iterdir())+list((R/'02_模型输入').iterdir())","files=list(RAW.iterdir())"),
    ("R/'02_模型输入/", "RAW/'"),
    ("R/'01_原始输入/", "RAW/'")])
emit(diet/'02_diet_core_pretest.py','00b_diet_pca.py',[
    ("sys.path.insert(0,'E:/new article/26_核心关联复核与投稿修订/_runtime')",'# Dependencies are supplied by the project environment.'),
    ("O=Path(__file__).parent;P=Path('E:/new article/37_审稿意见实质修订与复审/新增分析复现归档/输入/03_分析与结果/association_models')",PATHS+"\nP=A;O=D"),
    ("P/'single_genus_all_associations.tsv'","FROZEN/'single_genus_all_associations.tsv'")])

sex=ARCH/'04_性别与高读数复核/复现代码'
emit(sex/'02_sex_contrasts.py','09_sex_contrasts.py',[
    ('from common import RAW,X,XD,Z,Y,TT,BG,BC,glm,effect,contrast,bh,scale,np,pd,sm,norm','from common import RAW,X,XD,Z,Y,TT,BG,BC,glm,effect,contrast,bh,scale,np,pd,sm,norm,ARCH,RUN'),
    ("R=Path(__file__).resolve().parents[1];O=R/'新增分析结果'","R=RUN;O=RUN/'sex_results';O.mkdir(parents=True,exist_ok=True)"),
    ("R/'数据与代码/原R03_性别背景交互.tsv'","RUN/'results/R03_性别背景交互.tsv'")])
emit(sex/'08_abundance_influence.py','10_abundance_influence.py',[
    ('from common import X,XD,Z,Y,TT,BG,BC,np,pd,sm,bh','from common import X,XD,Z,Y,TT,BG,BC,np,pd,sm,bh,RUN'),
    ("R=Path(__file__).resolve().parents[1];O=R/'新增分析结果'","R=RUN;O=RUN/'sex_results';O.mkdir(parents=True,exist_ok=True);(R/'修订记录').mkdir(exist_ok=True)"),
    ("depth=pd.read_csv(R/'数据与代码/原39_named_depth.tsv',sep='\\t',index_col=0,dtype=str)","depth=pd.read_pickle(RUN/'results/_all_depth_family.pkl')")])

(ROOT/'reports').mkdir(exist_ok=True)
(ROOT/'reports/deepening_portability_changes.json').write_text(json.dumps(CHANGES,ensure_ascii=False,indent=2),encoding='utf-8')
print(f'Ported {len(CHANGES)} scripts without modifying archived sources.')
