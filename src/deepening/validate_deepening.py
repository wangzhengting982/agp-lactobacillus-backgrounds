"""Compare newly fitted outputs with submitted workbook and archived controls."""
from pathlib import Path
import os,json,math,sys
import pandas as pd,numpy as np,openpyxl
ROOT=Path(__file__).resolve().parents[2]
RUN=Path(os.environ.get('AGP_DEEPENING_OUTPUT',ROOT/'runs/deepening'))
ARCH=ROOT/'data/archive'
workbook=openpyxl.load_workbook(ROOT/'reference/submitted_results.xlsx',read_only=True,data_only=True)
checks=[];details=[]

def compare(name,new,old,source):
    info=dict(table=name,reference=source,new_shape=list(new.shape),reference_shape=list(old.shape))
    if new.shape!=old.shape or list(new.columns)!=list(old.columns):
        info.update(passed=False,reason='shape_or_column_order_difference');checks.append(info);return
    errors=[];max_abs=0.;max_rel=0.;numeric_cells=0
    for c in old:
        a=new[c];b=old[c]
        if pd.api.types.is_numeric_dtype(a) and not pd.api.types.is_bool_dtype(a):
            aa=pd.to_numeric(a,errors='coerce').to_numpy(float);bb=pd.to_numeric(b,errors='coerce').to_numpy(float)
            # Iterative fitting differences are permitted only below tight
            # numerical tolerance. Tiny nonzero p values require relative agreement.
            atol=1e-8;rtol=1e-7
            good=np.isclose(aa,bb,rtol=rtol,atol=atol,equal_nan=True)
            if c.startswith(('p','q_')) and c not in ['parameters','people','positive_people']:
                tiny=np.isfinite(bb)&(np.abs(bb)>0)&(np.abs(bb)<1e-7)
                good[tiny]=np.isclose(aa[tiny],bb[tiny],rtol=1e-6,atol=1e-300)
            finite=np.isfinite(aa)&np.isfinite(bb)
            if finite.any():
                diff=np.abs(aa[finite]-bb[finite]);max_abs=max(max_abs,float(diff.max()))
                max_rel=max(max_rel,float(np.max(diff/np.maximum(np.abs(bb[finite]),1e-300))))
            numeric_cells+=int(finite.sum())
        else:
            def canon(x):
                if x is None or (isinstance(x,float) and np.isnan(x)):return '<missing>'
                if isinstance(x,(bool,np.bool_)):return str(bool(x))
                return str(x)
            good=np.array([canon(x)==canon(y) for x,y in zip(a,b)])
        for i in np.flatnonzero(~good):
            errors.append(dict(row=int(i+2),column=str(c),new=str(a.iloc[i]),reference=str(b.iloc[i])))
    info.update(passed=not errors,numeric_cells=numeric_cells,max_absolute_error=max_abs,max_relative_error=max_rel,mismatch_count=len(errors))
    checks.append(info)
    if errors:details.append(dict(table=name,mismatches=errors[:50]))

for sheet in workbook.sheetnames:
    if sheet.split('_')[0] not in {'A01','A02','J01','J02','J03','L01','P01','P02','R01','R02','R03','R04','R05','S01','S02','S03','V01','V02','V03','H01','H02','H03','M01'}:continue
    values=list(workbook[sheet].values);old=pd.DataFrame(values[1:],columns=values[0])
    if sheet.startswith('M01'):
        meta=json.loads((RUN/'sex_results/M01_额外协变量可用性.json').read_text(encoding='utf8'));rows=[]
        # The submitted M01 contains these five fields; collection timestamps
        # were excluded when the supplement was condensed (15 audit-only rows).
        for field,entry in meta.items():
            if field not in ['contraceptive','pregnant','sex','age_years','diet_type']:continue
            for val,count in entry['counts'].items():
                by=entry.get('by_sex') or {};rows.append({'字段':field,'原值':val,'人数':count,'女性':by.get(val,{}).get('female'),'男性':by.get(val,{}).get('male')})
        new=pd.DataFrame(rows)
        new.to_csv(RUN/'sex_results/M01_额外协变量可用性.tsv',sep='\t',index=False)
    else:
        directory='sex_results' if sheet.startswith(('A02','H')) else 'results'
        path=RUN/directory/(sheet+'.tsv')
        if not path.exists():checks.append(dict(table=sheet,passed=False,reason='new_output_missing'));continue
        new=pd.read_csv(path,sep='\t')
    compare(sheet,new,old,'reference/submitted_results.xlsx')

# Verify reconstructed PCA scores/loads and food/nutrient mappings independently.
for label in ['food_PC','nutrient_PC']:
    for suffix in ['_得分_内部复现.tsv','_全部载荷.tsv']:
        name=label+suffix
        compare(name,pd.read_csv(RUN/'diet_preprocessing'/name,sep='\t'),pd.read_csv(ARCH/'02_统计输入/饮食预处理'/name,sep='\t'),'archived preprocessing control')
for name in ['_food_primary.pkl','_nutrient_primary.pkl']:
    new=pd.read_pickle(RUN/'diet_preprocessing'/name);old=pd.read_pickle(ARCH/'02_统计输入/饮食预处理'/name)
    checks.append(dict(table=name,passed=new.equals(old),new_shape=list(new.shape),reference_shape=list(old.shape),reference='archived intermediate control'))

report=dict(all_passed=all(c['passed'] for c in checks),checks=checks,mismatch_details=details,
    scope='All listed result tables are produced by new model fits; PCA starts from food and nutrient input tables. Baseline design, CLR and target counts come from runs/preparation freshly reconstructed inputs. R06/R07 taxonomy assignments and longitudinal index remain frozen inputs; taxonomy classification is not rerun here.',
    tolerance='absolute 1e-8 and relative 1e-7 for ordinary numeric cells; p/q below 1e-7 require relative 1e-6; strings, missingness, dimensions and row/column order exact')
(ROOT/'reports/deepening_validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ['checks','mismatch_details']},ensure_ascii=False))
for c in checks:print(c['table'],c['passed'],c.get('mismatch_count',''),c.get('max_absolute_error',''))
if details:print(json.dumps(details,ensure_ascii=False,indent=2)[:10000])
sys.exit(0 if report['all_passed'] else 1)
